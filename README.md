# README [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18868509.svg)](https://doi.org/10.5281/zenodo.18868509)

## Background
Quantum mechaincal phase space is important in semiclassical quantum mechanics, as it possesed a well-defined classical limit; however, exact sampling to the phase space function for anharmonic potentials has been a persistent challenge for decades with the Fourier transform in the definition. The current approximation schemes are, to varying degrees, restricted to Wigner-Kirkwood expansions, or to local harmonic approximations. 

## What we did 
In this work, we developed an *ab initio* method to build the thermal equilibrium Wigner $\mathcal{P}$-function and Husimi $\mathcal{Q}$-function using Ceperley's classical open polymer model. We also solved the equation of motion for quantum phase point with the correct deformation at thermal equilibrium. To see the mathematics and physics in this project, please refer to the .pdf document in this repo. 

## Reproduction
```bash
vi src/inputs/parameters.json # modify the parameters 
cd src/auto/
chmod +x submit_job.sh 
./submit_job.sh # run the simulation and visualize the result
