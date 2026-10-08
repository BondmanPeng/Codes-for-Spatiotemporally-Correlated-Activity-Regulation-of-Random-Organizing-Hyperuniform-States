# Codes for Spatiotemporally Correlated Activity Regulation of Random-Organizing Hyperuniform States

This repository contains C++ simulation codes and source data for **“Spatiotemporally correlated activity regulation of random-organizing hyperuniform states.”**

The simulations investigate how spatial and temporal correlations in driver-particle dynamics affect the density fluctuations of particles whose activity is locally regulated. Two types of drivers are considered: biased random-organization (BRO) particles and ideal-gas (IG) particles.

## Repository Structure

```text
source/
└── regulation_activity/
    ├── source_codes/
    │   ├── BRO/                 # Simulations with BRO drivers
    │   └── IG/                  # Simulations with ideal-gas drivers
    ├── available_data/
    │   ├── BRO/                 # Structure factors and plateau-fit data
    │   └── IG/
    │       ├── sq_SA_vary_DA/
    │       ├── sq_SBB_delta_combined/
    │       └── sq_SBB_minus_h0_T_h/
    └── README.txt               # Detailed dataset descriptions
```

## Physical Model

The system consists of two particle species:

- **Species A: driver particles.** Their dynamics generate the fluctuating local density field.
- **Species B: regulated particles.** Their activity depends on the local density of species A within a circular sensing region.

The two simulation models differ in the dynamics of species A. The BRO model uses biased random organization, while the IG model uses ideal-gas motion.

The supplied analysis codes and datasets characterize species-resolved structure factors, particle motion, and regulated-particle activity. The structure-factor data include the excess contribution

\[
\Delta S_B(q)=S_B(q;h)-S_B(q;h=0),
\]

where \(h\) controls the strength of activity regulation.

## Simulation and Analysis Codes

| Directory | File | Purpose |
| --- | --- | --- |
| `source_codes/BRO/` | `main_binary_randA.cpp` | Simulation with BRO drivers |
| `source_codes/BRO/` | `compute_sq_species.cpp` | Species-resolved structure factors |
| `source_codes/IG/` | `main_binary_idealA.cpp` | Simulation with ideal-gas drivers |
| `source_codes/IG/` | `compute_sq_species.cpp` | Species-resolved structure factors |
| `source_codes/IG/` | `compute_msd_A.cpp` | Mean-squared displacement of species A |
| `source_codes/IG/` | `compute_msd_B.cpp` | Mean-squared displacement of species B |
| `source_codes/IG/` | `compute_activity_B.cpp` | Activity of species B |

Supporting header files define particle properties, simulation dynamics, and numerical utilities.

## Main Parameters

| Parameter | Description |
| --- | --- |
| `rho_A`, `rho_B` | Densities of the driver and regulated particles |
| `h` | Strength of activity regulation |
| `h0` | Baseline activity scale |
| `R_A` | Radius of the circular sensing region |
| `c` | Correlation parameter for BRO driver kicks |
| `T` | IG simulation input related to driver diffusivity by \(D_A=T/2\) |
| `steps` | Number of simulation steps |
| `n_snap` | Number of requested snapshots |

## Source Data

The `available_data/` directory contains CSV tables for selected structure-factor curves, excess fluctuations, and plateau fits.

- **BRO data:** structure factors and excess fluctuations for sweeps over `c` and `h`, together with fitted plateau values and fit summaries.
- **IG data:** driver structure factors for different \(D_A\), regulated-particle structure factors, excess fluctuations, and associated fit results.

See [`source/regulation_activity/README.txt`](source/regulation_activity/README.txt) for detailed dataset definitions and parameter values.

## Code Requirements

The source code uses C++17. The BRO structure-factor analysis requires OpenMP; the IG CMake configuration detects OpenMP optionally.

The supplied CMake configurations reference test files that are not included in this package. Building through CMake therefore requires supplying those files or disabling the corresponding test targets. Plotting scripts referenced in the per-model README files are also not included.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
