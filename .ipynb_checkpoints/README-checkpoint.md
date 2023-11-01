## Abstract

The human gut microbiota is a complex ecosystem that affects a range of human physiology.  In order to explore the dynamics of the human gut microbiota, we used a system of ordinary differential equations to model mathematically the biomass of three microorganism populations: _Bacteroides thetaiotaomicron_, _Eubacterium rectale_, and _Methanobrevibacter smithii_. Additionally, we modeled the concentrations of relevant nutrients necessary to sustain these populations over time. Our model highlights the interactions and the competition among these three species. These three microorganisms were specifically chosen due to the system's end product, butyrate, which is a short chain fatty acid that aids in developing and maintaining the intestinal barrier in the human gut. The basis of our mathematical model assumes the gut is structured such that bacteria and nutrients exit the gut at a rate proportional to its volume, the rate of volumetric flow, and the biomass or concentration of the particular population or nutrient. We performed global sensitivity analyses using Sobol' sensitivities to estimate the relative importance of model parameters on simulation results.

## Repository Structure

Our main results can be found in ```model_julia.ipynb```.

/images: folder containing the png images that appear in the manuscript

This repository contains code files for the paper https://arxiv.org/abs/2303.12026

## Dependencies 

In this project, we utilized Julia v1.8.0. We list our package dependencies below, which are included in the Project.toml file in this repository.

[c7e460c6] ArgParse v1.1.4
[336ed68f] CSV v0.10.4
[13f3f980] CairoMakie v0.8.13
[1130ab10] DiffEqParamEstim v1.26.0
[0c46a032] DifferentialEquations v7.4.0
[31c24e10] Distributions v0.25.75
[af5da776] GlobalSensitivity v2.1.1
[4138dd39] JLD v0.13.2
[033835bb] JLD2 v0.4.23
[b964fa9f] LaTeXStrings v1.3.0
[429524aa] Optim v1.7.3
[7f7a1694] Optimization v3.9.0
[36348300] OptimizationOptimJL v0.1.2
[1dea7af3] OrdinaryDiffEq v6.27.2
[91a5bcdd] Plots v1.33.0
[8a4e6c94] QuasiMonteCarlo v0.2.9
[ed01d8cd] Sobol v1.5.0
[fce5fe82] Turing v0.21.12
[ade2ca70] Dates
[8bb1440f] DelimitedFiles

## BibTeX Citation

@misc{adrian2023,
      title={A mathematical model of B. thetaiotaomicron, M. smithii, and E. rectale interactions in the human gut}, 
      author={Melissa A. Adrian and Bruce P. Ayati and Ashutosh K. Mangalam},
      year={2023},
      eprint={2303.12026},
      archivePrefix={arXiv},
      primaryClass={q-bio.PE}
}
