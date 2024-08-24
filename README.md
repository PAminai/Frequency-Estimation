# ADSP: Frequency Estimation and CRB

## Author: Pouya Aminaie  


## Overview

This repository contains the MATLAB code and documentation for the first homework assignment in Advanced Digital Signal Processing (ADSP). The primary objective of this assignment is to estimate the frequency of a signal in the presence of noise using different techniques and analyze the performance of with different level of SNR using the Cramer-Rao Bound (CRB) and Mean Squared Error (MSE).

## Contents

- **Section 1: Theory**
  - Discusses the theoretical background required for noise generation and frequency estimation.
  - Details the computation of the sample covariance matrix, generation of random vectors, and the use of Power Spectral Density (PSD) for frequency estimation.
  
- **Section 2: Solutions**
  - **2.1 Noise Generation**
    - Explanation of noise generation using MATLAB’s `mvrnd` function.
    - Procedure for calculating the covariance matrix and MSE.
    - Discussion of the Monte Carlo method used in the simulations.
  - **2.2 Frequency Estimation**
    - Steps for estimating the frequency of a single sinusoidal signal.
    - Calculation of the covariance matrix in different scenarios.
    - Explanation of MATLAB functions and their roles in the implementation.
    - Results are presented showing MSE vs. SNR and CRB vs. SNR for various parameters.
  - **2.3 Frequency Modulation**
    - Analysis of the frequency modulation process and its effect on PSD.
  - **2.4 Arbitrary Modulation**
    - Strategy for detecting harmonics using PSD, with examples provided.

- **Section 3: Conclusion**
  - Summary of the work, key results, and possible future work, such as estimating spatial frequency in radar systems.

- **Appendix**
  - Contains figures and additional explanations that summarize the processes of noise generation and frequency estimation.

## Requirements

- MATLAB (tested with version R2023a)
- Signal Processing Toolbox

## How to Use

Download the following codes and reports from this repository:

1. [ADSP_report.pdf](https://github.com/PAminai/Frequency-Estimation/blob/main/ADSP_report.pdf)
2. [MyCode_HW1.m](https://github.com/PAminai/Frequency-Estimation/blob/main/MyCode_HW1.m)
3. [MyCode_HW1_1p3.m](https://github.com/PAminai/Frequency-Estimation/blob/main/MyCode_HW1_1p3.m)
4. [MyCode_HW1_1p4.m](https://github.com/PAminai/Frequency-Estimation/blob/main/MyCode_HW1_1p4.m)

