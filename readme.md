# Integrated Sensing, Communication, and Computing for Energy Efficient RIS-aided Wireless Federated Learning

Simulation code for the paper: **"Integrated Sensing, Communication, and Computing for Energy Efficient RIS-aided Wireless Federated Learning"**, which has been published in IEEE Transactions on Vehicular Technology.

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
  volume={75},
  number={3},
  pages={4336--4351},
  month={Mar.},
  year={2026},
  publisher={IEEE}
}
```

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

---

## 🤝 Contact

- Pouya Hosseini: hosseini.pouya7279@gmail.com

---
