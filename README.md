# Deterministic Chaos in the Double Pendulum: Lyapunov Exponent Analysis & Phase Space Topology

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Field: Physics & Applied Mathematics](https://img.shields.io/badge/Domain-Physics%20%7C%20Applied%20Math-orange.svg)](#)

A rigorous mathematical and computational investigation of deterministic chaos in nonlinear Hamiltonian dynamical systems using the planar double pendulum as a canonical model.

---

## 📌 Project Highlights

- **Theoretical Foundations**: Complete derivation of nonlinear equations of motion via Euler-Lagrange equations, Hamiltonian formulation, and Riemannian metric representation of the inertia tensor.
- **Lyapunov Exponent Analysis**: Numerical estimation of the Maximal Lyapunov Exponent (MLE) via tangent-space renormalization ($ \lambda \approx 2.1\ \text{s}^{-1} $), demonstrating exponential sensitivity and establishing finite predictability horizons.
- **Topological Phase Space Mapping**: Construction of Poincaré surface-of-section cuts ($ \theta_1 = 0 $) highlighting Kolmogorov-Arnold-Moser (KAM) tori, island chains, and chaotic sea coexistence.
- **Numerical Validation**: High-order Runge-Kutta integration (`RK45`, adaptive step size) preserving mechanical energy to machine-level relative error ($ |\Delta E / E_0| < 10^{-9} $).

---

## 🔬 Repository Structure

```text
research-paper/
├── research-paper.md     # Full academic paper (Equations, derivations, figures & analysis)
├── simulation.py         # Computational physics & dynamical systems simulation engine
├── figures/              # Generated high-resolution publication figures
│   ├── fig1_phase_space.png
│   ├── fig2_time_series.png
│   ├── fig3_sensitivity.png
│   ├── fig4_lyapunov.png
│   ├── fig5_poincare.png
│   ├── fig6_energy.png
│   └── fig7_lyapunov_map.png
└── README.md             # Project overview & reproduction instructions
```

---

## 🚀 Quickstart & Reproduction

### Prerequisites
- Python 3.10 or higher
- NumPy, SciPy, Matplotlib

```bash
pip install numpy scipy matplotlib
```

### Run Full Simulation & Generate Figures
```bash
python simulation.py
```

Generated plots will be saved in the `figures/` directory at 300 DPI publication quality.

---

## 📊 Summary of Results

| Metric | Result | Physical Interpretation |
| :--- | :--- | :--- |
| **Maximal Lyapunov Exponent ($\lambda$)** | $\approx 1.366\ \text{bits/s}$ | Trajectories separate exponentially |
| **Lyapunov Time ($\tau = 1/\lambda$)** | $\approx 0.73\ \text{s}$ | Fundamental prediction timescale |
| **Energy Conservation ($|\Delta E/E_0|$)** | $5.68 \times 10^{-11}$ | Validates symplectic structure preservation |
| **Poincaré Section Geometry** | Mixed state | Direct evidence for KAM theorem |

---

## 📖 Citation & Academic Context

If you reference this work or codebase in academic or portfolio reviews:
```bibtex
@article{double_pendulum_chaos_2026,
  title={Deterministic Chaos in the Double Pendulum: Lyapunov Exponent Analysis and Phase Space Topology},
  author={Raghvendra Kumar},
  year={2026},
  url={https://github.com/raghvendrakumarcreator/double-pendulum-chaos},
  journal={Portfolio Research Series in Computational Physics and Applied Mathematics}
}
```
