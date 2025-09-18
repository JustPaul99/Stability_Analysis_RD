

# 📊 Spatiotemporal Stability Analysis in MATLAB

## 🧠 Purpose

This project provides a MATLAB implementation for analyzing the **stability of spatiotemporal systems** using constrained least squares fitting and spectral analysis. It is designed to work with datasets that evolve over time and space, specifically reaction-diffusion systems, to help determine the stability of different Turing patterns. 

We provide the code for how we simulated the data, the code for how we analysed the data and the codes that we used for the experiments upon multiple datasets under different settings.

Codes created and based on the work....

---

## 📦 Repository Contents

```text
📁 Stability_Analysis_RD/Codes
 ┣ 📄 analyze_data.m                       # Core analysis function for extracting stability metrics from simulation data
 ┣ 📄 simulate_Klausmeier.m                # Euler-based simulation of the Klausmeier model (with noise and spatial dynamics)
 ┣ 📄 simulate_Klausmeier_RK.m             # RK4-based simulation of the Klausmeier model (higher accuracy time integration)
 ┣ 📄 simulate_and_analyze_combined.m      # Main script that runs simulations and analyzes results; generates raw data and basic plots
 ┗ 📄 plot_datasets_neat.m                 # Publication-grade plotting script using saved results; produces annotated figures and LaTeX-ready summaries
📄 README.md                               # Project documentation and workflow guide
```

---
## 🔍 Function: `Simulate_Klausmeier`

## 🧬 Data Generation Methods used for the Extended-Klausmeier Model

This repository supports two distinct approaches for simulating spatiotemporal data from the Klausmeier model:

### 1. **Euler Integration**  
Implemented in `simulate_Klausmeier.m`, this method uses a forward Euler scheme to evolve the system over time. It is simple and efficient for exploring general dynamics, especially when high precision is not critical. This method is more common and excepted for working with SDE's, as it represents an Euler-mayurama scheme.

### 2. **Runge-Kutta Integration (RK4)**  
Implemented in `simulate_Klausmeier_RK.m`, this method uses a fourth-order Runge-Kutta scheme for more accurate time integration. It is recommended when studying fine-scale dynamics or when numerical stability is important. With regards to implementation of SDEs, this is a bit more unconventional.

Both methods support different types of noise (`White`, `Correlated`, or `Uniform`) and allow flexible parameterization for exploring various ecological or physical regimes.

You can choose either method depending on your accuracy needs and computational budget. Both output matrices `u` and `v` that are compatible with the `analyze_data` function for downstream stability analysis.

Absolutely! Here's a concise and clear section you can add to your `README.md` to describe the **inputs and outputs** of your data simulation functions (`simulate_Klausmeier` and `simulate_Klausmeier_RK`):

markdown
---

## 📥 Inputs & 📤 Outputs

### Simulation Functions

Both `simulate_Klausmeier.m` and `simulate_Klausmeier_RK.m` generate spatiotemporal data based on the Klausmeier model. They accept the following inputs:

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
| `NoiseScale`         | Scaling factor for noise amplitude |

### Outputs

Both functions return:

- `u`: Matrix of vegetation density over time and space
- `v`: Matrix of water concentration over time and space

These outputs are compatible with the `analyze_data` function for downstream stability analysis.


---
## 🔍 Function: `analyze_data`

### 📌 Description

The function `analyze_data` takes in spatiotemporal data matrices `u` and `v`, extracts derivatives, fits a linear model under constraints, and performs a stability analysis by evaluating eigenvalues of a system matrix across spatial frequencies.

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
| `conf`             | Confidence intervals for each parameter |
| `StabilityMatrix`  | Function handle for the system matrix as a function of wavenumber `k` |
| `eigenvalue_function` | Function handle for the dominant eigenvalue |
| `k_max`            | Wavenumber that maximizes the dominant eigenvalue |
| `max_eigenvalue`   | Maximum eigenvalue at `k_max` (used to assess stability) |

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
- The dominant eigenvalue is computed across a range of `k` values.
- If `max_eigenvalue > 0`, the system is unstable; otherwise, it is stable.

---

## 🧪 Example Usage

```matlab
% Load or generate your data
[u, v] = data

% Define domain and sampling parameters
t = 10; L = 5;
tbegin = 1; tend = 100; tstepsize = 2;
xbegin = 1; xend = 50; stepsize = 1;

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

This repository includes two scripts for generating and analyzing data from the Klausmeier model. For the variations we used in the different experiments slight adaptations are neccessary of the code to make it functional, but if needed, can be sent upon request.

### 1. **Simulation & Analysis Script**
**File:** `simulate_and_analyze_combined.m`  
This script combines the simulation (`simulate_Klausmeier.m`) and analysis (`analyze_data.m`) functions to:
- Run multiple simulations across parameter ranges
- Extract fitted stability parameters from each run
- Save raw data and basic plots for inspection

It produces `.mat` files containing vegetation and water matrices (`u`, `v`), fitted parameters (`theta_data`), and eigenvalue curves (`plot_data`) for downstream use. It also plots the first version of the plots presented in the paper.

### 2. **Publication Plotting Script**
**File:** `plot_datasets_neat.m`  
This script loads the saved results from the first script and generates:
- Statistical summaries of eigenvalue curves
- Density overlays and confidence ellipses
- Annotated figures matching those presented in the paper

It is designed to produce clean, interpretable visualizations with the extensive focus on the stability parameters found for the seperate datasets. These plots should be similar to the ones presented in the paper.

---

## 📬 Contact

For questions or collaboration, reach out via GitHub Issues or email p.a.sanders@uu.nl.
```
