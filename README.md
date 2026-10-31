
# Numerical Simulation and Chaos Analysis of the Lorenz System

## Overview

This project investigates the Lorenz system, a classical nonlinear dynamical system known for exhibiting chaotic behaviour and sensitivity to initial conditions.

Using Python-based numerical simulations and mathematical analysis, this project explores the Lorenz attractor, equilibrium points, stability, Lyapunov exponents, and numerical bifurcation behaviour.

## Mathematical Model

The Lorenz system is described by the following differential equations:

dx/dt = σ(y - x)

dy/dt = x(ρ - z) - y

dz/dt = xy - βz

The classical parameter values used in this project are:

- σ = 10
- ρ = 28
- β = 8/3

## Project Objectives

1. Simulate and visualize the Lorenz attractor.
2. Investigate sensitivity to initial conditions.
3. Study the effect of parameter variation.
4. Calculate and visualize equilibrium points.
5. Analyse equilibrium stability using the Jacobian matrix and eigenvalues.
6. Estimate the largest Lyapunov exponent.
7. Study the variation of the Lyapunov exponent with respect to ρ.
8. Generate a numerical bifurcation diagram.

## Project Structure

```text
Lorenz_System_Simulation/
│
│Readme.txt
│   ├── 01_lorenzSystemSimulation.py
│   ├── 02_equilibrium_analysis.py
│   ├── 03_equilibrium_visualization.py
│   ├── 04_stability_analysis.py
│   ├── 05_lyapunov_analysis.py
│   ├── 06_lyapunov_convergence.py
│   ├── 07_lyapunov_vs_rho.py
│   └── 08_bifurcation_analysis.py
│
└── results/
    ├── figures/
    └── data/
```

## Experiments and Analysis

### 1. Lorenz System Simulation
Numerical integration of the Lorenz equations and visualization of the characteristic three-dimensional Lorenz attractor.

### 2. Sensitivity to Initial Conditions
Comparison of trajectories with slightly different initial conditions to demonstrate sensitive dependence on initial states.

### 3. Parameter Variation
Investigation of how changes in system parameters affect the resulting trajectories and dynamical behaviour.

### 4. Equilibrium Point Analysis
Calculation of equilibrium points of the Lorenz system and their corresponding mathematical properties.

### 5. Stability Analysis
Construction of the Jacobian matrix and analysis of its eigenvalues to determine the stability of equilibrium points.

### 6. Lyapunov Exponent Analysis
Numerical estimation of the largest Lyapunov exponent to investigate the divergence of nearby trajectories.

### 7. Lyapunov Exponent vs. Rho
Study of how the estimated largest Lyapunov exponent changes with the parameter ρ.

### 8. Bifurcation Analysis
Generation of a numerical bifurcation-style diagram using local maxima of the z variable for different values of ρ.

## Technologies Used

- Python
- NumPy
- SciPy
- Matplotlib

## Results

The project generates numerical results and visualizations related to:

- Three-dimensional Lorenz attractor
- Sensitivity to initial conditions
- Parameter variation
- Equilibrium points
- Stability of equilibrium points
- Lyapunov exponent estimates
- Lyapunov exponent variation with ρ
- Numerical bifurcation diagram

Generated figures are stored in `results/figures/`, while numerical CSV files are stored in `results/data/`.

## Mathematical Notes

For ρ > 1, the Lorenz system has two nonzero equilibrium points:

E± = (±√(β(ρ - 1)), ±√(β(ρ - 1)), ρ - 1)

The Jacobian matrix is:

J = [ -σ       σ       0  ]
    [ ρ-z     -1      -x  ]
    [  y       x      -β  ]

For the classical parameter regime, the nonzero equilibrium points lose stability at the Hopf threshold:

ρH = σ(σ + β + 3)/(σ - β - 1)

This stability threshold should not be interpreted as the exact onset of chaos.

## Limitations

The numerical results depend on the integration method, time step, initial conditions, transient removal, and simulation duration.

Lyapunov exponent values are numerical estimates and may vary with computational settings.

The bifurcation diagram is a numerical visualization based on local maxima and is not, by itself, a rigorous classification of all dynamical regimes.

## Authors

- Akriti kumari
- Rupam Kumari

Undergraduate Engineering Project

## Acknowledgement

This project was undertaken as an academic study of nonlinear dynamical systems and deterministic chaos.


