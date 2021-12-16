# SparseEcon

This code is Copyright Andreas Schaab and Allen T. Zhang, 2021-, and it is made available under the MIT license.

## Citation

If you find this software helpful, please cite our work:
- Schaab, A. and A. T. Zhang. Dynamic Programming in Continuous Time with Adaptive Sparse Grids. Working Paper.

## Overview and Applications
This repository provides a toolbox to solve dynamic programming problems in continuous time using adaptive sparse grids. Our method applies to a wide range of dynamic programming applications across many fields in economics.

We illustrate the power and applicability of our method by providing simple and accessible code to solve many dynamic programming applications using adaptive sparse grids. These can be found under ``local_repo_path/use_cases/``. Enumerated examples include:
- 01: A variant of Huggett (1993) with earnings risk modeled either as a two-state Markov chain or as a diffusion process
- 02: A two-asset variant of Aiyagari (1994), with a liquid asset (bond) and an illiquid asset (capital)
- 03: A two-asset portfolio choice problem with a riskfree asset (bond) and a risky asset (stocks)
- 04: Transition dynamics in response to MIT shocks in variants of Aiyagari (1994)
- 05: A global solution of Krusell-Smith (1998)
- 06: A growth model with non-convex production technologies in the spirit of Skiba (1978)
- 07: Life-cycle / OLG models, including variants with one-asset and two-asset portfolio choice problems
- 08: coming soon
- 09: A one-asset Heterogeneous Agent New Keynesian model in the spirit of McKay-Nakamura-Steinsson (2016)
- 10: A two-asset Heterogeneous Agent New Keynesian model in the spirit of Kaplan-Moll-Violante (2018)

In ``local_repo_path/use_cases/00_paper_replications``, we include code that replicates all output from both our working paper "_Dynamic Programming in Continuous Time with Adaptive Sparse Grids_" and the Handbook Chapter we wrote together with Johannes Brumm, Christopher Krause and Simon Scheidegger, "_Sparse grids for dynamic economic models_".

## Installation and Usage
Clone the repository and add the local path to your Matlab code, i.e.,
````
addpath(genpath('local_repo_path/lib/'))
````

## Requirements
- The software currently runs on MATLAB_R2021a

## Unit Tests
To run unit tests and verify reproducibility of use cases, run the following shell commands:
````
§ export PATH=$PATH:/Applications/MATLAB_R2021a.app/bin
§ cd local_repo_path
§ zsh run_tests.command
````
(Cross-platform test make file to come.)

## Acknowledgements
- We are deeply grateful to Ben Moll. We have learned and benefited greatly from his work and, in particular, the code he has made available online (https://benjaminmoll.com/codes/). If you find our software useful, please also consider citing: Achdou, Y., J. Han, J.-M. Lasry, P.-L. Lions, and B. Moll. Income and wealth distribution in macroeconomics: A continuous-time approach. Review of Economic Studies, Forthcoming. 2021.
- Some of our applications feature Linear Complementarity Problems, for which we use an LCP solver written by (Copyright, 2008) Yuval Tassa.
- All other references can be found in our Working Paper. 
- Also check out the GitHub repository for our Handbook Chapter (with J. Brumm, C. Krause and S. Scheidegger): https://github.com/SparseGridsForDynamicEcon.
