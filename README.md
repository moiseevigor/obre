# OBRE: Optimal Bias Robust Estimator Algorithms

[![DOI](https://zenodo.org/badge/5762/moiseevigor/obre.svg)](https://zenodo.org/badge/latestdoi/5762/moiseevigor/obre)

## Overview

This repository contains MATLAB and Fortran implementations of **Optimal Bias Robust Estimator (OBRE)** algorithms for extreme value analysis. The software is designed for robust parameter estimation of the Generalized Pareto Distribution (GPD), which is fundamental in modeling extreme events such as extreme wind speeds, floods, financial losses, and other rare but impactful phenomena.

## Theoretical Background

### Extreme Value Theory

Extreme Value Theory (EVT) provides the mathematical foundation for analyzing rare events that occur at the extremes of probability distributions. The theory is built on asymptotic results that describe the limiting behavior of extreme order statistics.

#### The Generalized Pareto Distribution

The **Generalized Pareto Distribution (GPD)** is central to the Peaks Over Threshold (POT) method in extreme value analysis. According to the **Pickands-Balkema-de Haan theorem**, for a wide class of underlying distributions, the conditional distribution of excesses over a high threshold converges to a GPD.

The GPD is defined as:

```
F(x) = 1 - [1 + ξ(x-μ)/σ]^(-1/ξ)    for ξ ≠ 0
F(x) = 1 - exp[-(x-μ)/σ]              for ξ = 0
```

where:
- **μ**: Location parameter (threshold)
- **σ > 0**: Scale parameter (α in this implementation)
- **ξ**: Shape parameter (β in this implementation)

### Parameter Estimation Methods

This implementation provides three parameter estimation approaches:

#### 1. Method of Moments (MoM)

The Method of Moments estimates parameters by equating sample moments to theoretical moments:

```
α̂ = (μ̂²/σ̂² + 1) × μ̂/2
β̂ = (μ̂²/σ̂² - 1)/2
```

where μ̂ and σ̂² are the sample mean and variance of excesses.

#### 2. Probability Weighted Moments (PWM)

PWM estimators, introduced by Hosking and Wallis (1987), are based on probability weighted moments:

```
M₀ = E[X]
M₁ = E[X × F(X)]
```

The PWM estimators are:

```
α̂ = 2 × M₀ × M₁ / (M₀ - 2M₁)
β̂ = M₀/(M₀ - 2M₁) - 2
```

#### 3. Optimal Bias Robust Estimators (OBRE)

OBRE methodology, developed by Ronchetti (1982, 1987), provides robust parameter estimation that is less sensitive to outliers while maintaining high efficiency under the assumed model.

**Ronchetti's Algorithm** consists of five steps:

1. **Initialization**: Start with PWM estimates
2. **Fisher Information**: Compute the Fisher information matrix J
3. **Eigendecomposition**: Find eigenvalues and eigenvectors of J
4. **Iterative Refinement**: Update parameter estimates using robust score functions
5. **Convergence**: Iterate until convergence criteria are met

The OBRE approach minimizes the maximum asymptotic bias over a neighborhood of the assumed model, making it robust against model misspecification and outliers.

## Applications

### Wind Engineering

The primary application domain is **structural wind engineering**, where accurate estimation of extreme wind speeds is crucial for:

- Building design and safety assessment
- Determination of design wind loads
- Estimation of return period wind speeds (VREF)

The return period formula implemented is:

```
VREF = μ + (α/β) × [1 - (n/(M×T))^(-β)]
```

where:
- **μ**: Threshold wind speed
- **n**: Number of exceedances
- **M**: Number of years of record
- **T**: Return period (years)

### Other Applications

- **Hydrology**: Flood frequency analysis
- **Finance**: Risk management and Value-at-Risk estimation
- **Insurance**: Catastrophic loss modeling
- **Seismology**: Earthquake magnitude analysis

## Software Implementation

### MATLAB Implementation

The MATLAB version provides a complete implementation with the following key functions:

#### Main Functions

- **`obre.m`**: Main function implementing the complete OBRE algorithm
- **`meth1.m`**: Method of Moments parameter estimation
- **`meth2.m`**: Probability Weighted Moments parameter estimation
- **`step345.m`**: Steps 3-5 of Ronchetti's algorithm

#### Supporting Functions

- **`fishinfo.m`**: Fisher information matrix calculation
- **`jacobi.m`**: Eigenvalue decomposition
- **`goodoffit.m`**: Goodness-of-fit testing (Anderson-Darling, Cramér-von Mises)
- **`quangp.m`**: Quantile estimation for GPD
- **`weight.m`**: Robust weight function calculation
- **`covmat.m`**: Asymptotic covariance matrix estimation

#### Utility Functions

- **`readinfo.m`**: Data preprocessing and structure initialization
- **`vmean.m`**, **`vstdev.m`**: Statistical calculations
- **`dgsint.m`**: Numerical integration
- **`score1.m`**, **`score2.m`**: Score function calculations

### Fortran Implementation

The Fortran version provides high-performance implementation suitable for intensive computational applications (see `fortran/` directory).

## Usage

### Basic Usage

```matlab
% Basic OBRE analysis with default parameters
[alpha, beta, nex, threshold, vref_obre, vref_mom, vref_pwm] = obre(25.8, 50);
```

### Advanced Usage

```matlab
% OBRE analysis with custom parameters
config = struct(...
    'incr', 0.2, ...      % Threshold increment
    'rep', 5, ...         % Number of threshold repetitions
    'c', 1.5, ...         % Bound on influence function
    'epsiln', 0.001, ...  % Convergence tolerance
    'factor', 10.0 ...    % Integration speed factor
);

[alpha, beta, nex, threshold, vref_obre, vref_mom, vref_pwm, weights, errors] = ...
    obre(25.8, 50, config);
```

### Input Parameters

- **`initialthreshold`**: Initial threshold value for analysis
- **`returnperiod`**: Return period for VREF calculation (years)
- **`DATINIT`** (optional): Configuration structure with fields:
  - `incr`: Threshold increment for sensitivity analysis
  - `rep`: Number of threshold repetitions
  - `c`: Bound on the influence function (robustness parameter)
  - `epsiln`: Convergence tolerance for iterative algorithm
  - `factor`: Integration factor (1.0 = slow, 50.0 = fast)
  - `sim`: Number of simulations for goodness-of-fit testing

### Output Parameters

- **`alpha`**: Scale parameter estimates (vector)
- **`beta`**: Shape parameter estimates (vector)
- **`nex`**: Number of exceedances for each threshold
- **`threshold`**: Final threshold values used
- **`vref_obre`**: VREF estimates using OBRE
- **`vref_mom`**: VREF estimates using Method of Moments
- **`vref_pwm`**: VREF estimates using PWM
- **`weights`**: Robust weights assigned to observations
- **`errors`**: Structure containing standard errors and confidence intervals

## Data Requirements

The software expects MATLAB `.mat` files containing:
- **Peak values**: Vector of observed extreme values
- **Independence**: Data should be processed for temporal independence
- **Stationarity**: Underlying process should be stationary

### Data Preprocessing

For wind speed analysis, data preprocessing typically includes:
1. **Peak extraction**: Identify local maxima above threshold
2. **Independence**: Ensure temporal independence between peaks
3. **Stationarity**: Verify stationarity of the underlying process

## Mathematical Details

### Robustness Properties

The OBRE estimators satisfy:

1. **Breakdown Point**: 1/n (same as maximum likelihood)
2. **Bounded Influence**: Influence function is bounded by parameter `c`
3. **High Efficiency**: Asymptotic efficiency close to maximum likelihood under the model

### Goodness-of-Fit Testing

The implementation provides several goodness-of-fit tests:

#### Anderson-Darling Test
```
A² = -n - (1/n) × Σᵢ[(2i-1) × ln(uᵢ) + (2n+1-2i) × ln(1-uᵢ)]
```

#### Cramér-von Mises Test
```
W² = (1/12n) + Σᵢ[uᵢ - (2i-1)/(2n)]²
```

where uᵢ are the probability integral transforms.

## Literature References

### Core References

1. **Ronchetti, E.** (1982). *Robust testing in linear models*. PhD Thesis, ETH Zurich.

2. **Ronchetti, E.** (1987). *Robustness aspects of model choice*. Statistica Sinica, 7, 327-338.

3. **Pickands, J.** (1975). *Statistical inference using extreme order statistics*. Annals of Statistics, 3, 119-131.

4. **Balkema, A.A. and de Haan, L.** (1974). *Residual life time at great age*. Annals of Probability, 2, 792-804.

### Methodology References

5. **Hosking, J.R.M. and Wallis, J.R.** (1987). *Parameter and quantile estimation for the generalized Pareto distribution*. Technometrics, 29, 339-349.

6. **Davison, A.C. and Smith, R.L.** (1990). *Models for exceedances over high thresholds*. Journal of the Royal Statistical Society, Series B, 52, 393-442.

7. **Coles, S.** (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer-Verlag.

### Applications References

8. **Simiu, E. and Heckert, N.A.** (1996). *Extreme wind distribution tails: a peaks over threshold approach*. Journal of Structural Engineering, 122, 539-547.

9. **Cook, N.J.** (1985). *The Designer's Guide to Wind Loading of Building Structures*. Butterworths.

### Robust Statistics References

10. **Huber, P.J. and Ronchetti, E.** (2009). *Robust Statistics*, 2nd Edition. John Wiley & Sons.

11. **Maronna, R.A., Martin, R.D., and Yohai, V.J.** (2006). *Robust Statistics: Theory and Methods*. John Wiley & Sons.

## Software Information

### Development History

- **Original Fortran Implementation**: Joanna E. Mills (TUNS, 1996), modified by D.J. Dupuis (1997) and Andrew P. Wilson (1999)
- **MATLAB Port**: Igor Moiseev (2011)
- **Version**: 0.1 for MATLAB

### License

This software is distributed under the **GNU General Public License v3.0**:

```
OBRE for Matlab
Copyright (C) 2011, Igor Moiseev

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
```

### System Requirements

- **MATLAB**: R2008a or later
- **Fortran**: Standard Fortran 90/95 compiler
- **Memory**: Depends on dataset size (typically < 1GB for most applications)
- **Storage**: Minimal (algorithms are computationally intensive, not storage intensive)

## Troubleshooting

### Common Issues

1. **Convergence Problems**: 
   - Reduce `epsiln` tolerance
   - Increase `factor` for faster integration
   - Check data for extreme outliers

2. **Negative Scale Parameters**:
   - Indicates poor fit or inappropriate threshold
   - Try different threshold values
   - Verify data preprocessing

3. **High Shape Parameter (β > 1)**:
   - May indicate heavy-tailed distribution
   - Check for data quality issues
   - Consider alternative distributions

### Performance Optimization

- Use `factor` parameter to balance accuracy vs. speed
- For large datasets, consider parallel processing
- Pre-process data to remove obvious outliers

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

### Development Priorities

- GPU acceleration for large datasets
- Additional robust estimators (M-estimators, S-estimators)
- Automated threshold selection methods
- Extended goodness-of-fit testing suite

## Citation

If you use this software in your research, please cite:

```bibtex
@misc{moiseev2011obre,
  title={OBRE: Optimal Bias Robust Estimator Algorithms},
  author={Moiseev, Igor},
  year={2011},
  doi={10.5281/zenodo.xxxxx},
  url={https://github.com/moiseevigor/obre}
}
```

## Contact

For questions, bug reports, or collaboration inquiries, please open an issue on GitHub or contact the maintainers.

---

*This software implements state-of-the-art robust statistical methods for extreme value analysis, providing reliable parameter estimation in the presence of model uncertainty and outliers.*
