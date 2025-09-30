

# 📊 Spatiotemporal Stability Analysis in MATLAB

## 🧠 Purpose

This project provides a MATLAB implementation for analyzing the **stability of spatiotemporal systems** using constrained least squares fitting and spectral analysis to create dispersion relations. It is designed to work with datasets that evolve over time and space, and is developed specifically for reaction-diffusion systems. 

We provide the code for how we simulated the data, the code for how we analysed the data and the code that we used for the statistical analysis of the multiple experimental setups.

Codes created for the work....

---

## 📦 Repository Contents

```text
📁 Stability_Analysis_RD/Codes
 ┣ 📄 analyze_data.m                       # Core analysis function for extracting stability metrics from simulation data
 ┣ 📄 simulate_Klausmeier.m                # Euler-based simulation of the Klausmeier model (with noise and spatial dynamics)
 ┣ 📄 simulate_and_analyze_combined.m      # Main script that runs simulations and analyzes results; generates raw data and basic plots
 ┗ 📄 plot_datasets_neat.m                 # plotting script to get similar plots to the paper.
📄 README.md                               # Project documentation and workflow guide
```

---
## 🔍 Function: `Simulate_Klausmeier`

### 📌 Description

`simulate_Klausmeier.m` generates spatiotemporal data based on the Extended Klausmeier model as defined in the paper. This method uses an Euler-Mayurama scheme with three different types of noise (`White`, `Correlated`, or `Uniform`).

### 📎 Syntax
```matlab
[u, v] = simulate_Klausmeier(delta, h, m, pvec, t, tsteps, L, xsteps, noiseType, correlation_length, NoiseScale)
```
### 📥 inputs

| Parameter            | Description |
|----------------------|-------------|
| `delta`              | Diffusion coefficient for vegetation (`v`) |
| `h`                  | Strength of negative feedback in vegetation growth |
| `m`                  | Mortality rate of vegetation |
| `pvec`               | Time-dependent precipitation vector |
| `t`                  | Total simulation time |
| `tsteps`             | Number of time steps |
| `L`                  | Length of the spatial domain |
| `xsteps`             | Number of spatial steps |
| `noiseType`          | Type of noise: `'White'`, `'Correlated'`, or `'Uniform'` |
| `correlation_length` | Spatial correlation length for noise |
| `NoiseScale`         | Scaling factor for noise strength |

### 📤 Outputs
- `u`: Matrix of vegetation density over time and space
- `v`: Matrix of water concentration over time and space

These outputs are compatible with `analyze_data`.


---
## 🔍 Function: `analyze_data`

### 📌 Description

The function `analyze_data` takes in spatiotemporal data matrices `u` and `v`, extracts derivatives, fits a linear model under constraints, and performs a stability analysis by evaluating eigenvalues of a system matrix across different spatial eigenmodes.

### 📎 Syntax

```matlab
[theta, conf, StabilityMatrix, eigenvalue_function, k_max, max_eigenvalue] = analyze_data(u, v, t, L, tbegin, tend, tstepsize, xbegin, xend, stepsize)
```

### 📥 Inputs

| Parameter      | Description |
|----------------|-------------|
| `u`, `v`       | 2D matrices of spatiotemporal data (e.g., velocity fields) |
| `t`            | Total time duration of the dataset |
| `L`            | Length of the spatial domain |
| `tbegin`, `tend` | Start and end indices for time sampling |
| `tstepsize`    | Time sampling interval |
| `xbegin`, `xend` | Start and end indices for spatial sampling |
| `stepsize`     | Spatial sampling interval |

### 📤 Outputs

| Output             | Description |
|--------------------|-------------|
| `theta`            | Fitted model parameters |
| `conf`             | Gaussian estimated confidence intervals for each parameter |
| `StabilityMatrix`  | Function handle for the system matrix as a function of wavenumber `k` |
| `eigenvalue_function` | Function handle for the eigenvalue for different eigenmodes |
| `k_max`            | Dominant wavenumber for the found model |
| `max_eigenvalue`   | Eigenvalue at the dominant wavenumber |

---

## ⚙️ How It Works

### 1. **Data Sampling**
The function samples a subset of the full dataset based on user-defined time and space windows.

### 2. **Derivative Computation**
- Time derivatives are computed using forward differences.
- Spatial second derivatives are computed using central differences.

### 3. **Model Fitting**
- A linear model is fit to the data using `lsqcurvefit`.
- Constraints include bounds, equality conditions, and inequality conditions to reflect physical assumptions.

### 4. **Stability Analysis**
- A system matrix is constructed as a function of spatial frequency `k`.
- The dominant eigenmode `k_max` and its corresponding eigenvalue `max_eigenvalue` are calculated from this matrix.

---

## 🧪 Example Usage

```matlab
%%optional data creation
%Domain setup
t = 1; L = 40;
tsteps=10000;
xsteps=400;

% Setup parameters for klausmeier
delta=0.01; % Turing
h=0.1;
m=0.5;
pbegin = 6; 
pend = 6;
pvec = pbegin + (pend - pbegin) * (0:dt:t) / t;

%Noise properties
correlation_length=0.1;
noiseType='Correlated';
NoiseScale=1; %Noise Strength

% Load or generate your data
[u, v] = simulate_Klausmeier(delta, h, m, pvec, t, tsteps, L, xsteps, noiseType, correlation_length, NoiseScale)
%[u,v] = Your_data

% Define domain and sampling parameters for analysis
tbegin = 0; tend = tsteps; tstepsize = 1;
xbegin = 0; xend = xsteps; stepsize = 1;

% Run analysis
[theta, conf, StabilityMatrix, eigenvalue_function, k_max, max_eigenvalue] = ...
    analyze_data(u, v, t, L, tbegin, tend, tstepsize, xbegin, xend, stepsize);

% Optional: Print stability result
if max_eigenvalue > 0
    fprintf('Unstable: λ = %.2f at k = %.2f\n', max_eigenvalue, k_max);
else
    fprintf('Stable: λ = %.2f at k = %.2f\n', max_eigenvalue, k_max);
end
```
---

## 📈 Workflow Overview

This repository includes two scripts for creating multiple datasets from the Extended Klausmeier model that are then being analyzed for stability. One specific experimental setup has been provided for the creating of the datasets for different noise situations. For other experimental setups, you can follow the values in Table 1 of the paper and make different loops for these settings.

### 1. **Simulation & Analysis Script**
**File:** `simulate_and_analyze_combined.m`  
This script combines the simulation (`simulate_Klausmeier.m`) and analysis (`analyze_data.m`) functions to:
- Run multiple simulations across different experimentsettings
- Extract fitted stability parameters from each run
- Save raw data and basic plots for inspection

It produces `.mat` files containing vegetation and water matrices (`u`, `v`), fitted parameters (`theta_data`), and eigenvalue curves (`plot_data`) for downstream use. It also plots the first version of the plots presented in the paper.

### 2. **Publication Plotting Script**
**File:** `plot_datasets_neat.m`  
This script loads the saved results from the previous script and plots it in a similar way to the original paper.

---

## 📬 Contact

For questions or collaboration, reach out via GitHub Issues or email p.a.sanders@uu.nl.
```
