# gut_microbiota

## Abstract

The human gut microbiota is a complex ecosystem that affects a range of human physiology.  In order to explore the dynamics of the human gut microbiota, we used a system of ordinary differential equations to model mathematically the biomass of three microorganism populations: _Bacteroides thetaiotaomicron_, _Eubacterium rectale_, and _Methanobrevibacter smithii_. Additionally, we modeled the concentrations of relevant nutrients necessary to sustain these populations over time. Our model highlights the interactions and the competition among these three species. These three microorganisms were specifically chosen due to the system's end product, butyrate, which is a short chain fatty acid that aids in developing and maintaining the intestinal barrier in the human gut. The basis of our mathematical model assumes the gut is structured such that bacteria and nutrients exit the gut at a rate proportional to its volume, the rate of volumetric flow, and the biomass or concentration of the particular population or nutrient. We performed global sensitivity analyses using Sobol' sensitivities to estimate the relative importance of model parameters on simulation results.

## Repository Structure

Our main results can be found in ```model_julia.ipynb```.

/images: folder containing the png images that appear in the manuscript

This repository contains code files for the paper https://arxiv.org/abs/2303.12026

https://drive.google.com/drive/folders/1DIVYCq_e2JhSJtS0IqshmzGxZoACegO3

## BibTeX Citation

@misc{adrian2023,
      title={A mathematical model of B. thetaiotaomicron, M. smithii, and E. rectale interactions in the human gut}, 
      author={Melissa A. Adrian and Bruce P. Ayati and Ashutosh K. Mangalam},
      year={2023},
      eprint={2303.12026},
      archivePrefix={arXiv},
      primaryClass={q-bio.PE}
}
