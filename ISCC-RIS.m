%--------------------------------------------------------------------------
% Clear command window, workspace variables, and close all figures
%--------------------------------------------------------------------------
clc;
clear;
close all;

%--------------------------------------------------------------------------
% Monte Carlo Simulation Settings
% If this variable is true, the simulation runs for num_iterations.
% Otherwise, it will run only once to show convergence.
%--------------------------------------------------------------------------
run_monte_carlo = false;

if run_monte_carlo
    num_iterations = 50;  % Number of Monte Carlo simulation iterations
else
    num_iterations = 1;   % Single run to observe convergence
end

%--------------------------------------------------------------------------
% Pre-allocating memory for storing results
% This improves the loop execution speed.
%--------------------------------------------------------------------------
results_EE = zeros(num_iterations, 1);    % Store Energy Efficiency results
results_Rate = zeros(num_iterations, 1);  % Store Total Rate results
results_Power = zeros(num_iterations, 1); % Store Total Power results

%--------------------------------------------------------------------------
% Start of the main simulation loop
%--------------------------------------------------------------------------
tic; % Start timer
for iter = 1:num_iterations
    fprintf('=========== Monte Carlo Iteration: %d ===========\n', iter);

    %----------------------------------------------------------------------
    % Define main system parameters
    %----------------------------------------------------------------------
    M = 16; % Number of Base Station (BS) antennas
    N = 40; % Number of Reconfigurable Intelligent Surface (RIS) elements
    K = 4;  % Number of users
    T = 3;  % Number of targets for sensing

    % Local Computing (Federated Learning) Parameters
    fk_min = 2.5e9; % Minimum user processing frequency (Hz)
    fk_max = 4.0e9; % Maximum user processing frequency (Hz)
    mu_k = 1e-28;   % Effective switched capacitance coefficient of the processor
    Qk = 20e6;      % Number of computation samples (Cycles)
    Ck = 20;        % Number of computation cycles per sample
    I = 8;          % Number of local computation iterations

    % Power and Noise Parameters
    sigma2 = 10^(-80/10) / 1000; % Additive White Gaussian Noise (AWGN) power in Watts
    P_C = 10^(25/10) / 1000;     % Constant power consumption of BS circuits in Watts
    P_th_loc = 10^(35/10) / 1000;% Maximum local power consumption threshold for a user in Watts
    P_max = 10^(40/10) / 1000;   % Maximum BS transmit power in Watts
    xi = 1.1;                    % Inverse of the energy conversion coefficient at the transmitter

    % Communication Parameters
    B = 20e6;       % Total system bandwidth (Hz)
    tau_k_req = 10^(3/10); % Minimum required Signal-to-Interference-plus-Noise Ratio (SINR) for a user

    % Sensing Parameters
    epsilon = 0.1; % Maximum allowed error between designed and ideal beampattern (MSE)
    theta_bar = [-35, 0, 35] * pi/180; % Central angles of the targets (radians)
    delta_theta = 10 * pi/180;       % Angular width of each target (radians)

    % Channel Model Parameters (Path Loss and Rician)
    C0_dB = -30;   % Path loss at the reference distance (dB)
    d_BR = 50;     % Distance between BS and RIS (meters)
    d_BU = 60;     % Distance between BS and user (meters)
    d_RU = 8;      % Distance between RIS and user (meters)
    alpha_BR = 2.5;% Path loss exponent for BS-RIS link
    alpha_BU = 3.5;% Path loss exponent for BS-User link
    alpha_RU = 2.5;% Path loss exponent for RIS-User link
    beta_BR = 10^(3/10); % Rician factor for BS-RIS link
    beta_RU = 10^(0/10); % Rician factor for RIS-User link
    beta_BU = 10^(3/10); % Rician factor for BS-User link

    % Calculate linear path loss from dB value
    L_BR = 10^(C0_dB/10) * (d_BR)^(-alpha_BR);
    L_BU = 10^(C0_dB/10) * (d_BU)^(-alpha_BU);
    L_RU = 10^(C0_dB/10) * (d_RU)^(-alpha_RU);

    %----------------------------------------------------------------------
    % Generate Communication Channels based on the Rician model
    %----------------------------------------------------------------------
    % Channel between BS and RIS
    % Initialize channel matrices
    G = zeros(N, M);
    H_BU = zeros(M, K);
    H_RU = zeros(N, K);

    % Parameters
    alpha = 0.5; % CSI quality

    % Generate channels
    for m = 1:N
        for n = 1:M
            % True G channel (scalar per element)
            G_true = sqrt(L_BR) * (sqrt(1 / (1 + beta_BR)) * (randn + 1i * randn) / sqrt(2));
            % Error term (scalar)
            E = (randn + 1i * randn) / sqrt(2);
            % Imperfect CSI for G
            G(m, n) = alpha * G_true + sqrt(1 - alpha^2) * E;
        end
    end

    for k = 1:K
        % True H_BU channel (M x 1)
        H_BU_true = sqrt(L_BU) * (sqrt(beta_BU / (1 + beta_BU)) * ones(M, 1) + ...
            sqrt(1 / (1 + beta_BU)) * (randn(M, 1) + 1i * randn(M, 1)) / sqrt(2));
        % Error term (M x 1)
        E_BU = (randn(M, 1) + 1i * randn(M, 1)) / sqrt(2);
        % Imperfect CSI for H_BU
        H_BU(:, k) = alpha * H_BU_true + sqrt(1 - alpha^2) * E_BU;

        % True H_RU channel (N x 1)
        H_RU_true = sqrt(L_RU) * (sqrt(1 / (1 + beta_RU)) * (randn(N, 1) + 1i * randn(N, 1)) / sqrt(2));
        % Error term (N x 1)
        E_RU = (randn(N, 1) + 1i * randn(N, 1)) / sqrt(2);
        % Imperfect CSI for H_RU
        H_RU(:, k) = alpha * H_RU_true + sqrt(1 - alpha^2) * E_RU;
    end

    %----------------------------------------------------------------------
    % Define the ideal beampattern for sensing
    %----------------------------------------------------------------------
    theta_grid = linspace(-pi/2, pi/2, 181); % Grid of angles for pattern calculation
    P_id = compute_ideal_beampattern(theta_grid, theta_bar, delta_theta);

    %----------------------------------------------------------------------
    % Initialize optimization variables
    %----------------------------------------------------------------------
    % W_init: Initial beamforming matrix (includes communication and sensing parts)
    W_init = sqrt(P_max / (2 * M * (K + M))) * (randn(M, K + M) + 1i * randn(M, K + M));
    % phi_init: Initial phases of RIS elements
    phi_init = 2 * pi * rand(N, 1);
    % f_init: Initial processing frequencies of users
    f_init = fk_min * ones(K, 1);
    % b_init: Initial bandwidth allocated to users
    b_init = B / K * ones(K, 1);
    % alpha_init: Initial user selection variables
    alpha_init = rand(K, 1);
    % gamma_init: Initial scaling factor for the beampattern
    gamma_init = 1.0;

    % Aggregate all optimization variables into a single vector x_k
    x_k = [real(W_init(:)); imag(W_init(:)); phi_init; f_init; b_init; alpha_init; gamma_init];
    
    % Define lower and upper bounds for the optimization variables
    num_W_real_vars = M * (K+M) * 2;
    lb_x = [-inf(num_W_real_vars, 1); zeros(N, 1); fk_min * ones(K, 1); zeros(K, 1); zeros(K, 1); 0];
    ub_x = [ inf(num_W_real_vars, 1); 2*pi*ones(N, 1); fk_max * ones(K, 1); B * ones(K, 1);  ones(K, 1); inf];

    %----------------------------------------------------------------------
    % Parameters for the SLP (Successive Linear Programming) Algorithm
    % These parameters control the trust region and penalty factor.
    %----------------------------------------------------------------------
    rho_bad = 0.25;    % Lower threshold for the ratio of actual to predicted reduction
    rho_good = 0.75;   % Upper threshold for the reduction ratio
    zeta = 2;          % Trust region radius update factor
    kappa = 2;         % Penalty factor increase factor
    epsilon1 = 1e-4;   % Tolerance for stopping condition (small step norm)
    Delta_k = 0.1;     % Initial trust region radius
    Delta_max = 1;     % Maximum trust region radius
    Delta_min = 1e-6;  % Minimum trust region radius
    lambda_p = 1e5;    % Penalty factor for binary variables alpha
    lambda_k = 1e10;   % Initial penalty factor for constraint violations
    lambda_max = 1e12; % Maximum penalty factor

    %----------------------------------------------------------------------
    % Start of the main SLP algorithm loop
    %----------------------------------------------------------------------
    max_slp_iter = 150; % Maximum number of SLP iterations
    ee_history = zeros(max_slp_iter, 1); % To store the EE convergence history

    for r = 1:max_slp_iter
        % 1. Unpack variables from the vector x_k
        [Wv, Psiv_diag, fv, bv, alphav, gammav] = unpack_variables(x_k, M, K, N);
        
        % 2. Build linearized constraints using first-order Taylor expansion
        [A_ineq, b_ineq] = build_linearized_constraints(Wv, Psiv_diag, fv, bv, alphav, gammav, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K);
        
        % 3. Build the linearized objective function's gradient
        f_grad = build_linearized_objective_gradient(Wv, Psiv_diag, fv, bv, alphav, lambda_p, H_RU, H_BU, G, sigma2, mu_k, xi, P_C, M, N, K);
        
        num_vars = length(x_k);
        num_ineq_constraints = size(A_ineq, 1);
        
        % 4. Define the bounds for the step d_k based on the trust region
        d_lb = max(-Delta_k, (lb_x - x_k));
        d_ub = min(Delta_k, (ub_x - x_k));
        
        % 5. Solve the Linear Programming (LP) subproblem using CVX
        % This subproblem finds the optimal step d_k within the trust region.
        % The objective is to minimize the sum of the linearized objective and the penalty for constraint violations.
        cvx_begin quiet
            variable d_k(num_vars) % Optimization step
            variable t(num_ineq_constraints) nonnegative % Slack variables for penalty

            minimize( f_grad' * d_k + lambda_k * sum(t) )
            subject to
                A_ineq * d_k - b_ineq <= t; % Linearized constraints
                d_k >= d_lb;                % Lower bound on step (trust region)
                d_k <= d_ub;                % Upper bound on step (trust region)
        cvx_end

        fval_lp = cvx_optval; % Optimal value of the subproblem's objective function
        
        % Check the status of the CVX solver
        if ~contains(cvx_status, 'Solved')
            fprintf('  SLP Iteration %d: Error - CVX could not solve the problem. Status: %s\n', r, cvx_status);
            % If the problem is not solved, shrink the trust region and continue
            Delta_k = max(Delta_k / (zeta*2), Delta_min);
            if Delta_k == Delta_min
                 fprintf('  SLP Iteration %d: Cannot recover. Stopping.\n', r);
                 break;
            end
            continue;
        end
        
        % Stopping condition: if the step size is very small, we are close to an optimum
        if norm(d_k) < epsilon1
            disp('Optimal solution found (norm of d_k is small).');
            break;
        end
        
        % 6. Calculate the Merit Function
        % This function combines the original objective function and the sum of constraint violations.
        [true_obj_k, const_viol_k, ee_k] = calculate_merit_function(x_k, lambda_p, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K, xi, P_C);
        ee_history(r) = ee_k; % Store EE value for plotting convergence
        
        x_new = x_k + d_k; % New candidate point
        [true_obj_new, const_viol_new] = calculate_merit_function(x_new, lambda_p, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K, xi, P_C);
        
        % 7. Update trust region and penalty factor
        actual_merit_k = true_obj_k + lambda_k * const_viol_k;
        actual_merit_new = true_obj_new + lambda_k * const_viol_new;
        
        actual_reduction = actual_merit_k - actual_merit_new; % Actual reduction in the merit function
        predicted_reduction = lambda_k * const_viol_k - fval_lp; % Predicted reduction by the linear model
        
        % Calculate the ratio rho to decide on the step
        if predicted_reduction < 1e-9
            rho = -1; % Avoid division by zero
        else
            rho = actual_reduction / predicted_reduction;
        end
        
        if ~run_monte_carlo
            fprintf('  SLP Iteration %d: |d|=%.3e, rho=%.3f, Delta=%.3e, EE=%.4f, Viol=%.3e\n', r, norm(d_k), rho, Delta_k, ee_k, const_viol_k);
        end
        
        % Update the point and trust region radius based on rho
        if rho > rho_bad % The step is good, accept it
            x_k = x_k + d_k;
            if rho > rho_good % The step is very good, expand the trust region
                Delta_k = min(zeta * Delta_k, Delta_max);
            end
        else % The step is bad, reject it and shrink the trust region
            Delta_k = max(Delta_k / zeta, Delta_min);
        end

        % Update the penalty factor lambda_k
        if const_viol_new > 1e-3 || rho < rho_bad
            lambda_k = min(kappa * lambda_k, lambda_max);
        end
    end
    
    %----------------------------------------------------------------------
    % Extract and store the final results
    %----------------------------------------------------------------------
    [W_final, Psi_final, f_final, b_final, alpha_final, gamma_final] = unpack_variables(x_k, M, K, N);
    
    [~, ~, final_EE, final_Rate, final_Power] = calculate_merit_function(x_k, lambda_p, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K, xi, P_C);
    
    results_EE(iter) = final_EE;
    results_Power(iter) = final_Power;
    results_Rate(iter) = final_Rate;
    
    fprintf('>>> Iteration %d Result: EE = %f (Rate=%.2f Mbps / Power=%.2f W)\n', iter, final_EE, final_Rate/1e6, final_Power);
end
toc; % End timer

%--------------------------------------------------------------------------
% Display Final Results
%--------------------------------------------------------------------------
if run_monte_carlo
    figure;
    boxplot(results_EE, 'Labels', {'Energy Efficiency (bit/J)'});
    title('Distribution of Energy Efficiency over Monte Carlo Runs');
    ylabel('Energy Efficiency (bit/J)');
    grid on;
else
    figure;
    plot(1:r, ee_history(1:r), '-o', 'LineWidth', 2);
    title('Energy Efficiency Convergence Plot');
    xlabel('SLP Iteration');
    ylabel('Energy Efficiency (bit/J)');
    grid on;
end

fprintf('\n\n========== Final Results (Average of %d runs) ==========\n', num_iterations);
fprintf('Average Energy Efficiency: %.4f bit/J\n', mean(results_EE, 'omitnan'));
fprintf('Average Normalized EE:     %.4f bit/J/Hz\n', mean(results_EE, 'omitnan')/B);
fprintf('Average Total Rate:        %.4f Mbps\n', mean(results_Rate, 'omitnan') / 1e6);
fprintf('Average Total Power:       %.4f Watts\n', mean(results_Power, 'omitnan'));
fprintf('==============================================================\n\n');

% Check if the final solution satisfies the constraints
check_final_constraints(x_k, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K, fk_min, fk_max);


%==========================================================================
% Helper Functions
%==========================================================================

function [W, Psi_diag, f, b, alpha, gamma] = unpack_variables(x, M, K, N)
    % This function converts the optimization vector x back to the original matrix and vector variables.
    
    % Extract the beamforming matrix W
    num_W_vars = M * (K+M);
    W_real_vec = x(1 : num_W_vars);
    W_imag_vec = x(num_W_vars + 1 : 2 * num_W_vars);
    W = reshape(W_real_vec + 1i * W_imag_vec, M, K+M);

    % Extract the diagonal matrix of RIS phases (Psi)
    phi_start = 2 * num_W_vars + 1;
    phi = x(phi_start : phi_start + N - 1);
    Psi_diag = diag(exp(1i * phi)); % Psi is a diagonal matrix with elements e^(j*phi_n)

    % Extract the processing frequencies vector f
    f_start = phi_start + N;
    f = x(f_start : f_start + K - 1);

    % Extract the bandwidth vector b
    b_start = f_start + K;
    b = x(b_start : b_start + K - 1);

    % Extract the user selection vector alpha
    alpha_start = b_start + K;
    alpha = x(alpha_start : alpha_start + K - 1);

    % Extract the scaling variable gamma
    gamma = x(end);
end

%--------------------------------------------------------------------------

function P_id = compute_ideal_beampattern(theta_grid, theta_bar, delta_theta)
    % This function constructs the ideal beampattern.
    % The pattern is 1 in the target directions and 0 elsewhere.
    P_id = zeros(size(theta_grid));
    for t = 1:length(theta_bar)
        P_id(theta_grid >= theta_bar(t) - delta_theta/2 & theta_grid <= theta_bar(t) + delta_theta/2) = 1;
    end
end

%--------------------------------------------------------------------------

function a = steering_vector(theta, M)
    % This function creates the steering vector for a given angle.
    delta_lambda = 0.5; % Antenna spacing in terms of half-wavelength
    m_idx = (0:M-1)';
    a = exp(1i * 2 * pi * delta_lambda * m_idx * sin(theta));
end

%--------------------------------------------------------------------------

function [SINR_val, All_grad_W, All_grad_phi] = compute_SINR_and_grads(W, Psi_diag, H_RU, H_BU, G, sigma2, K, M, N)
    % This function calculates the SINR values and their gradients with respect to W and phi.
    h_eff_all = H_RU' * Psi_diag * G + H_BU'; % Total effective channel
    
    SINR_val = zeros(K, 1);
    All_grad_W = cell(K, 1);
    All_grad_phi = zeros(N, K);
    
    for k = 1:K
        hk_eff = h_eff_all(k, :); % Effective channel for user k
        wk = W(:, k);             % Beamforming vector for user k
        
        % Calculate signal power
        signal_power = abs(hk_eff * wk)^2;
        
        % Calculate interference power
        interference_power = 0;
        all_indices = 1:(K+M);
        interferer_indices = all_indices(all_indices ~= k);
        for j = interferer_indices
            interference_power = interference_power + abs(hk_eff * W(:, j))^2;
        end
        total_interference = interference_power + sigma2;
        
        % Calculate SINR
        SINR_val(k) = signal_power / total_interference;
        
        % Calculate SINR gradient with respect to W
        grad_W_k = zeros(M, K+M, 'like', 1i);
        grad_num_Wk = 2 * hk_eff' * (hk_eff * wk);
        grad_W_k(:,k) = grad_num_Wk / total_interference;
        
        for j = interferer_indices
            grad_den_Wj = 2 * hk_eff' * (hk_eff * W(:, j));
            grad_W_k(:,j) = grad_W_k(:,j) - (signal_power / total_interference^2) * grad_den_Wj;
        end
        All_grad_W{k} = grad_W_k;
        
        % Calculate SINR gradient with respect to phi
        grad_phi_k = zeros(N, 1);
        psi_vec = diag(Psi_diag);
        for n=1:N
            d_hk_eff_n = H_RU(n,k)' * (1i * psi_vec(n)) * G(n,:); % Derivative of effective channel wrt phi_n
            d_num_phi_n = 2 * real( conj(hk_eff*wk) * (d_hk_eff_n * wk) );
            d_den_phi_n = 0;
            for j=1:(K+M)
                 d_den_phi_n = d_den_phi_n + 2*real( conj(hk_eff * W(:,j)) * (d_hk_eff_n * W(:,j)) );
            end
            grad_phi_k(n) = (d_num_phi_n * total_interference - signal_power * d_den_phi_n) / total_interference^2;
        end
        All_grad_phi(:,k) = grad_phi_k;
    end
end

%--------------------------------------------------------------------------

function f_obj = build_linearized_objective_gradient(Wv, Psiv_diag, fv, bv, alphav, lambda_p, H_RU, H_BU, G, sigma2, mu_k, xi, P_C, M, N, K)
    % This function builds the gradient of the objective function of problem P2 for the LP subproblem.
    % The objective function is a combination of negative energy efficiency and a penalty term for alpha.
    
    [SINRv, All_grad_Wv, All_grad_phiv] = compute_SINR_and_grads(Wv, Psiv_diag, H_RU, H_BU, G, sigma2, K, M, N);
    
    safe_SINRv = max(0, SINRv);
    log_term = log2(1 + safe_SINRv);

    % Calculate total rate (Rv) and total power (Pv) at the current point
    Rv = sum(alphav .* bv .* log_term);
    Pv = sum(alphav .* mu_k .* fv.^3) + xi * norm(Wv, 'fro')^2 + P_C;
    if abs(Pv) < 1e-9; Pv = 1e-9; end % Avoid division by zero

    % Calculate the gradient of Rate (R) with respect to different variables
    grad_R_W = zeros(M, K+M, 'like', 1i);
    grad_R_phi = zeros(N, 1);
    for k=1:K
       if (1+safe_SINRv(k)) > 1e-9
           common_term = (alphav(k) * bv(k)) / (log(2) * (1 + safe_SINRv(k)));
           grad_R_W = grad_R_W + common_term * All_grad_Wv{k};
           grad_R_phi = grad_R_phi + common_term * All_grad_phiv(:,k);
       end
    end
    grad_R_W_real = [real(grad_R_W(:)); imag(grad_R_W(:))];
    grad_R_f = zeros(K, 1);
    grad_R_b = alphav .* log_term;
    grad_R_alpha = bv .* log_term;
    
    % Calculate the gradient of Power (P) with respect to different variables
    grad_P_W_real = 2 * xi * [real(Wv(:)); imag(Wv(:))];
    grad_P_phi = zeros(N, 1);
    grad_P_f = alphav .* 3 .* mu_k .* fv.^2;
    grad_P_b = zeros(K, 1);
    grad_P_alpha = mu_k .* fv.^3;

    % Calculate the gradient of Energy Efficiency (EE = R/P) using the quotient rule
    grad_EE_W = (Pv * grad_R_W_real - Rv * grad_P_W_real) / Pv^2;
    grad_EE_phi = (Pv * grad_R_phi - Rv * grad_P_phi) / Pv^2;
    grad_EE_f = (Pv * grad_R_f - Rv * grad_P_f) / Pv^2;
    grad_EE_b = (Pv * grad_R_b - Rv * grad_P_b) / Pv^2;
    grad_EE_alpha = (Pv * grad_R_alpha - Rv * grad_P_alpha) / Pv^2;
    grad_EE_gamma = 0; % EE is not dependent on gamma
    
    % Gradient of the alpha penalty term
    grad_penalty_alpha = lambda_p * (1 - 2 * alphav);
    
    % Total gradient of the objective function (-EE + Penalty)
    grad_total_W = -grad_EE_W;
    grad_total_phi = -grad_EE_phi;
    grad_total_f = -grad_EE_f;
    grad_total_b = -grad_EE_b;
    grad_total_alpha = -grad_EE_alpha + grad_penalty_alpha;
    grad_total_gamma = -grad_EE_gamma;
    
    % Aggregate gradients into a single vector
    f_obj = [grad_total_W; grad_total_phi; grad_total_f; grad_total_b; grad_total_alpha; grad_total_gamma];
end

%--------------------------------------------------------------------------

function [A_ineq, b_ineq] = build_linearized_constraints(Wv, Psiv_diag, fv, bv, alphav, gammav, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K)
    % This function linearizes the constraints of problem P2 and puts them in the form A*d <= b.
    
    num_W_real_vars = M * (K+M) * 2;
    num_vars = num_W_real_vars + N + K + K + K + 1;
    
    % Define indices of variables in the total vector
    w_indices = 1:num_W_real_vars;
    phi_indices = num_W_real_vars+1 : num_W_real_vars+N;
    f_indices = num_W_real_vars + N + 1 : num_W_real_vars + N + K;
    b_indices = num_W_real_vars + N + K + 1 : num_W_real_vars + N + 2*K;
    alpha_indices = num_W_real_vars + N + 2*K + 1 : num_W_real_vars + N + 3*K;
    gamma_idx = num_vars;

    A_ineq = [];
    b_ineq = [];

    % Constraint C2: Max transmit power ||W||^2 <= P_max
    grad_c2_W = 2 * [real(Wv(:)); imag(Wv(:))];
    A_c2 = zeros(1, num_vars);
    A_c2(w_indices) = grad_c2_W';
    b_c2 = P_max - norm(Wv, 'fro')^2;
    A_ineq = [A_ineq; A_c2];
    b_ineq = [b_ineq; b_c2];

    % Constraint C5: Total allocated bandwidth sum(alpha_k * b_k) <= B
    A_c5 = zeros(1, num_vars);
    grad_c5_alpha = bv;
    grad_c5_b = alphav;
    A_c5(alpha_indices) = grad_c5_alpha';
    A_c5(b_indices) = grad_c5_b';
    b_c5 = B - sum(alphav .* bv);
    A_ineq = [A_ineq; A_c5];
    b_ineq = [b_ineq; b_c5];
    
    % Constraint C6: Number of selected users 1 <= sum(alpha_k) <= K
    grad_c6_lower = zeros(1, num_vars);
    grad_c6_lower(alpha_indices) = -ones(1, K);
    rhs_c6_lower = -1 + sum(alphav);

    grad_c6_upper = zeros(1, num_vars);
    grad_c6_upper(alpha_indices) = ones(1, K);
    rhs_c6_upper = K - sum(alphav);

    A_ineq = [A_ineq; grad_c6_lower; grad_c6_upper];
    b_ineq = [b_ineq; rhs_c6_lower; rhs_c6_upper];
    
    % Constraint C4: Local power consumption alpha_k * (mu_k * f_k^3 - P_th_loc) <= 0
    A_c4 = zeros(K, num_vars);
    b_c4 = zeros(K, 1);
    for k=1:K
        pk_loc_v = mu_k * fv(k)^3;
        g_k_v = pk_loc_v - P_th_loc;
        grad_g_k_f = 3 * mu_k * fv(k)^2;
        A_c4(k, alpha_indices(k)) = g_k_v;
        A_c4(k, f_indices(k)) = alphav(k) * grad_g_k_f;
        b_c4(k) = - (alphav(k) * g_k_v);
    end
    A_ineq = [A_ineq; A_c4];
    b_ineq = [b_ineq; b_c4];

    % Constraint C9: Minimum required SINR alpha_k * (tau_k_req - SINR_k) <= 0
    [SINRv, All_grad_Wv, All_grad_phiv] = compute_SINR_and_grads(Wv, Psiv_diag, H_RU, H_BU, G, sigma2, K, M, N);
    A_c9 = zeros(K, num_vars);
    b_c9 = zeros(K, 1);
    for k=1:K
        g_k_v = tau_k_req - SINRv(k);
        grad_g_W_k_complex = -All_grad_Wv{k};
        grad_g_W_k_real = [real(grad_g_W_k_complex(:)); imag(grad_g_W_k_complex(:))];
        grad_g_phi_k = -All_grad_phiv(:,k);
        
        A_c9(k, alpha_indices(k)) = g_k_v;
        A_c9(k, w_indices) = alphav(k) * grad_g_W_k_real';
        A_c9(k, phi_indices) = alphav(k) * grad_g_phi_k';
        b_c9(k) = - (alphav(k) * g_k_v);
    end
    A_ineq = [A_ineq; A_c9];
    b_ineq = [b_ineq; b_c9];
    
    % Constraint C8: Beampattern error E(gamma, W) <= epsilon
    [Ev, grad_E_W_real, grad_E_gamma] = compute_beampattern_error_and_grads(Wv, gammav, P_id, theta_grid, M);
    A_c8 = zeros(1, num_vars);
    A_c8(w_indices) = grad_E_W_real';
    A_c8(gamma_idx) = grad_E_gamma;
    b_c8 = epsilon - Ev;
    A_ineq = [A_ineq; A_c8];
    b_ineq = [b_ineq; b_c8];
end

%--------------------------------------------------------------------------

function [E, grad_W_real, grad_gamma] = compute_beampattern_error_and_grads(W, gamma, P_id, theta_grid, M)
    % This function calculates the Mean Squared Error (MSE) of the beampattern and its gradients.
    E = 0;
    grad_W = zeros(size(W), 'like', 1i);
    grad_gamma = 0;
    
    L = length(theta_grid);
    for l = 1:L
        a_l = steering_vector(theta_grid(l), M); % Steering vector at angle l
        designed_power = real(a_l' * W * W' * a_l); % Designed power
        ideal_power = gamma * P_id(l); % Scaled ideal power
        error_term = ideal_power - designed_power;
        E = E + error_term^2;
        
        % Calculate gradients
        grad_gamma = grad_gamma + 2 * error_term * P_id(l);
        grad_W = grad_W - 2 * error_term * (2 * a_l * (a_l' * W));
    end
    % Normalize by the number of angle grid points
    E = E / L;
    grad_gamma = grad_gamma / L;
    grad_W = grad_W / L;
    grad_W_real = [real(grad_W(:)); imag(grad_W(:))];
end

%--------------------------------------------------------------------------

function [obj, viol, ee_val, rate_val, power_val] = calculate_merit_function(x, lambda_p, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K, xi, P_C)
    % This function calculates the merit function value (objective + penalty) and the main performance metrics.
    [W, Psi_diag, f, b, alpha, gamma] = unpack_variables(x, M, K, N);
    
    % Calculate SINR for all users
    h_eff_all = H_RU' * Psi_diag * G + H_BU';
    SINR_vals = zeros(K,1);
    for k = 1:K
        hk_eff = h_eff_all(k, :);
        signal_power = abs(hk_eff * W(:, k))^2;
        all_indices = 1:(K+M);
        interferer_indices = all_indices(all_indices ~= k);
        interference_power = norm(hk_eff * W(:, interferer_indices))^2;
        SINR_vals(k) = signal_power / (interference_power + sigma2);
    end
    
    % Calculate main metrics: rate, power, and energy efficiency
    rate_val = sum(alpha .* b .* log2(1 + max(0, SINR_vals)));
    power_val = sum(alpha .* mu_k .* f.^3) + xi * norm(W, 'fro')^2 + P_C;
    if power_val < 1e-9; power_val = 1e-9; end
    ee_val = rate_val / power_val;
    
    % Calculate the main objective function (with penalty)
    penalty_alpha = lambda_p * sum(alpha .* (1 - alpha));
    obj = -ee_val + penalty_alpha;
    
    % Calculate the sum of constraint violations
    viol = 0;
    viol = viol + max(0, norm(W, 'fro')^2 - P_max);
    viol = viol + max(0, sum(alpha .* b) - B);
    viol = viol + max(0, 1 - sum(alpha));
    for k=1:K
        viol = viol + max(0, alpha(k) * (mu_k * f(k)^3 - P_th_loc));
        viol = viol + max(0, alpha(k) * (tau_k_req - SINR_vals(k)));
    end
    [E_val, ~, ~] = compute_beampattern_error_and_grads(W, gamma, P_id, theta_grid, M);
    viol = viol + max(0, E_val - epsilon);
end

%--------------------------------------------------------------------------

function check_final_constraints(x, H_RU, H_BU, G, sigma2, tau_k_req, mu_k, P_th_loc, P_max, B, epsilon, P_id, theta_grid, M, N, K, fk_min, fk_max)
    % This function checks whether the final solution satisfies the original problem constraints.
    [W, Psi_diag, f, b, alpha, gamma] = unpack_variables(x, M, K, N);
    tolerance = 1e-4; % A small tolerance for numerical comparisons

    fprintf('========== Final Constraint Check (Numbering based on P1 in the paper) ==========\n');
    
    c1_viol = any(f < fk_min - tolerance) || any(f > fk_max + tolerance);
    fprintf('C1: Frequency bounds        \t|\t %s\n', get_status(~c1_viol));

    c2_viol = norm(W, 'fro')^2 > P_max + tolerance;
    fprintf('C2: Max transmit power      \t|\t %s (Value: %.4f / Limit: %.4f)\n', get_status(~c2_viol), norm(W, 'fro')^2, P_max);

    c3_viol = any(abs(abs(diag(Psi_diag)) - 1) > tolerance);
    fprintf('C3: RIS phase shifts (unit-modulus) \t|\t %s\n', get_status(~c3_viol));

    local_powers = mu_k .* f.^3;
    c4_viol = false;
    for k_check = 1:K
        if alpha(k_check) > 0.5 && local_powers(k_check) > P_th_loc + tolerance
            c4_viol = true;
            fprintf('    (User %d violates C4: Power %.4f W > Limit %.4f W)\n', k_check, local_powers(k_check), P_th_loc);
            break;
        end
    end
    fprintf('C4: Local power consumption      \t|\t %s\n', get_status(~c4_viol));
    
    c5_viol = sum(alpha .* b) > B + tolerance;
    fprintf('C5: Total bandwidth         \t|\t %s (Value: %.4f / Limit: %.4f)\n', get_status(~c5_viol), sum(alpha .* b), B);

    c6_viol = (sum(alpha) < 1 - tolerance) || (sum(alpha) > K + tolerance);
    fprintf('C6: Number of selected users    \t|\t %s (Value: %.4f)\n', get_status(~c6_viol), sum(alpha));
    
    c7_viol = any(alpha < -tolerance) || any(alpha > 1 + tolerance);
    fprintf('C7: Alpha variable in [0,1]   \t|\t %s\n', get_status(~c7_viol));

    [E_val, ~, ~] = compute_beampattern_error_and_grads(W, gamma, P_id, theta_grid, M);
    c8_viol = E_val > epsilon + tolerance;
    fprintf('C8: Beampattern error          \t|\t %s (Value: %.4e / Limit: %.4e)\n', get_status(~c8_viol), E_val, epsilon);

    h_eff_all = H_RU' * Psi_diag * G + H_BU';
    SINR_vals = zeros(K,1);
    for k = 1:K
        hk_eff = h_eff_all(k, :);
        signal_power = abs(hk_eff * W(:, k))^2;
        all_indices = 1:(K+M);
        interferer_indices = all_indices(all_indices ~= k);
        interference_power = norm(hk_eff * W(:, interferer_indices))^2;
        SINR_vals(k) = signal_power / (interference_power + sigma2);
    end
    c9_viol = false;
    for k_check = 1:K
        if alpha(k_check) > 0.5 && SINR_vals(k_check) < tau_k_req - tolerance
            c9_viol = true;
            fprintf('    (User %d violates C9: SINR is %.4f < Required %.4f)\n', k_check, SINR_vals(k_check), tau_k_req);
            break;
        end
    end
    fprintf('C9: Minimum required SINR    \t|\t %s\n', get_status(~c9_viol));

    fprintf('=================================================================\n');
end

%--------------------------------------------------------------------------

function status_str = get_status(is_satisfied)
    % A simple function to display "OK" or "VIOLATED" status
    if is_satisfied
        status_str = 'OK';
    else
        status_str = '!!! VIOLATED !!!';
    end
end
