# ISCC-RIS: Integrated Sensing, Communication, and Computing for Energy Efficient RIS-aided Wireless Federated Learning

Simulation code for the paper: **"Integrated Sensing, Communication, and Computing for Energy Efficient RIS-aided Wireless Federated Learning"**, which has been submitted to IEEE Transactions on Vehicular Technology.

---

## 📄 About the Article

This work explores an energy-efficient integrated sensing, communication, and computing (ISCC) system for privacy-preserving wireless federated learning (FL) at edge devices in 6G and beyond wireless networks. We propose an efficient user selection method based on local computing time concerning high-power processing to enhance user participation and employ reconfigurable intelligent surface (RIS) to improve channel reliability in the wireless FL system.

To maximize the energy efficiency (EE), we leverage a nonlinear optimization problem as a mixed-integer nonlinear program (MINLP) with constraints on power, computing frequency, bandwidth, and target sensing metric based on beamforming design matrices. The MINLP is transformed into a tractable nonlinear program, approximated via linear problems using Taylor expansion and solved with an iterative algorithm based on successive linear programming (SLP).

**Key Contributions:**
- Advanced energy-efficient ISCC system leveraging RIS to improve SINR and optimize channel conditions
- Effective user selection method based on local computing time for low-latency users
- Joint optimization of beamforming matrix, RIS reflecting coefficients, local computing frequency, bandwidth allocation, and user selection
- MINLP formulation transformed into tractable linear problems using Taylor expansion
- SLP algorithm with trust region method and penalty factor for convergence guarantee

---

## 📝 Citation

If you use this code for research that results in publications, please cite our original article:

```bibtex
@article{kesargheh2025integrated,
  title={Integrated Sensing, Communication, and Computing for Energy Efficient RIS-aided Wireless Federated Learning},
  author={Kesargheh, Mohammad Mansour and Hosseini, Pouya and Nouri, Nima and Zahedi, Abdulhamid and Abouei, Jamshid and Mohammadi, Arash},
  journal={IEEE Transactions on Vehicular Technology},
  year={2025},
  note={Submitted}
}
```

---

## 📁 Code Structure

```
ISCC-RIS/
├── ISCC-RIS.m          # Main simulation script
└── README.md           # This file
```

### File Description

| File | Description |
|------|-------------|
| `ISCC-RIS.m` | Main MATLAB script containing the complete simulation framework, including channel generation, SLP optimization algorithm, helper functions, and result visualization |

---

## 🔧 Requirements

- **MATLAB** (R2020a or later recommended)
- **CVX** (Cvx Research, Inc.) - Required for solving the linear programming subproblems
  - Download from: http://cvxr.com/cvx/
  - Install and run `cvx_setup` before executing the code

---

## 🚀 Usage

### 1. Setup CVX

Ensure CVX is properly installed and initialized:

```matlab
cvx_setup
```

### 2. Run the Simulation

Simply execute the main script in MATLAB:

```matlab
ISCC-RIS
```

### 3. Simulation Modes

The script supports two simulation modes controlled by the `run_monte_carlo` variable:

| Mode | Setting | Description |
|------|---------|-------------|
| **Convergence** | `run_monte_carlo = false` | Single run showing EE convergence over SLP iterations |
| **Monte Carlo** | `run_monte_carlo = true` | Multiple runs (default: 50) for statistical analysis |

### 4. Key Parameters

Adjustable system parameters in the main script:

```matlab
M = 16;                    % Number of BS antennas
N = 40;                    % Number of RIS elements
K = 4;                     % Number of users
T = 3;                     % Number of sensing targets
B = 20e6;                  % Bandwidth (Hz)
P_max = 10^(40/10)/1000;   % Max transmit power (W)
sigma2 = 10^(-80/10)/1000; % Noise power (W)
```

---

## 📊 Output

### Single Run Mode (`run_monte_carlo = false`)

- **Figure 1**: Energy Efficiency convergence plot vs. SLP iterations
- **Console Output**: Iteration details including step norm, rho, Delta, EE, and constraint violations
- **Final Metrics**: Average EE, Rate, and Power consumption

### Monte Carlo Mode (`run_monte_carlo = true`)

- **Figure 1**: Boxplot showing EE distribution across all runs
- **Console Output**: Results for each iteration and final averages

### Constraint Check

The script automatically checks all constraints (C1-C9) against the original problem formulation and reports whether each is satisfied or violated.

---

## 🔍 Code Features

### Main Functions

| Function | Description |
|----------|-------------|
| `unpack_variables()` | Converts optimization vector back to matrix/vector variables |
| `compute_ideal_beampattern()` | Constructs ideal radiation pattern for target directions |
| `steering_vector()` | Creates steering vector for given angle |
| `compute_SINR_and_grads()` | Calculates SINR and gradients w.r.t. W and phi |
| `build_linearized_objective_gradient()` | Constructs gradient of linearized objective |
| `build_linearized_constraints()` | Linearizes constraints using first-order Taylor expansion |
| `compute_beampattern_error_and_grads()` | Computes beampattern MSE and gradients |
| `calculate_merit_function()` | Evaluates merit function with penalty |
| `check_final_constraints()` | Verifies all constraints on final solution |

### SLP Algorithm Parameters

```matlab
rho_bad = 0.25;      % Lower threshold for reduction ratio
rho_good = 0.75;     % Upper threshold for reduction ratio
zeta = 2;            % Trust region update factor
kappa = 2;           % Penalty factor increase factor
Delta_k = 0.1;       % Initial trust region radius
lambda_p = 1e5;      % Penalty for binary alpha variables
lambda_k = 1e10;     % Initial penalty factor
epsilon1 = 1e-4;     % Stopping tolerance
```

---

## 📈 System Model Overview

The proposed system consists of:

1. **BS** with M antennas serving K single-antenna users
2. **RIS** with N elements near users to enhance communication
3. **Sensing**: BS detects T targets with fixed LoS connections
4. **FL**: Users perform local computation and send updates to BS
5. **User Selection**: Based on local computing time criterion

### Optimization Variables

- Beamforming matrix W ∈ ℂ^(M×(K+M))
- RIS phase shifts Ψ = diag(ψ)
- Local computing frequency f_k
- Bandwidth allocation b_k
- User selection α_k
- Beampattern scaling factor γ

---

## ⚠️ Notes

1. **Channel Estimation**: The code includes both perfect and imperfect CSI cases with parameter α ∈ [0,1]
2. **Path Loss**: Rician fading model with specified path loss exponents and Rician factors
3. **Sensing**: Beampattern similarity metric with MSE constraint ε
4. **NOMA Support**: Code is structured for FDMA; NOMA extension is described in the paper

---

## 🤝 Contact

**Corresponding Authors:**
- Jamshid Abouei: abouei@yazd.ac.ir
- Arash Mohammadi: arash.mohammadi@concordia.ca

**First Author:**
- Mohammad Mansour Kesargheh: mohammadmansourkesargheh@gmail.com

---

## 📚 Related References

1. A. Adhikary, A. Deb Raha, Y. Qiao, W. Saad, Z. Han, and C. Seon Hong, "Holographic MIMO with integrated sensing and communication for energy-efficient cell-free 6G networks," *IEEE Internet of Things J.*, vol. 11, no. 19, pp. 30617-30635, Oct. 2024.
2. M. Rihan, A. Zappone, S. Buzzi, D. Wubben, and A. Dekorsy, "Energy efficiency maximization for active RIS-aided integrated sensing and communication," *J. Wireless Comm. Network*, vol. 2024, no. 20, pp. 1-22, Apr. 2024.
3. J. Nocedal and S. J. Wright, *Numerical Optimization*. Springer, Jul. 2006.

---

## 📜 License

This code is provided for academic and research purposes. Please cite the original paper if you use this code in your work.

---
