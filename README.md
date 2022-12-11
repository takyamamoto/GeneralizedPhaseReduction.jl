# GeneralizedPhaseReduction.jl

[![Build Status](https://travis-ci.com/takyamamoto/GeneralizedPhaseModel.jl.svg?branch=main)](https://travis-ci.com/takyamamoto/GeneralizedPhaseModel.jl)
[![Coverage](https://codecov.io/gh/takyamamoto/GeneralizedPhaseModel.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/takyamamoto/GeneralizedPhaseModel.jl)
[![Coverage](https://coveralls.io/repos/github/takyamamoto/GeneralizedPhaseModel.jl/badge.svg?branch=main)](https://coveralls.io/github/takyamamoto/GeneralizedPhaseModel.jl?branch=main)

Julia package for generalized phase reduction method.

<img src="https://raw.githubusercontent.com/takyamamoto/GeneralizedPhaseReduction.jl/master/figures/logo.png" width="30%" align="right" />


## Requirements
Julia above ver 1.6. 

### Installation
Run the following command. 
```julia
pkg > add https://github.com/takyamamoto/GeneralizedPhaseReduction.jl.git
```

### Dependent packages 
```
ForwardDiff.jl, LinearAlgebra.jl, Random.jl, Statistics.jl, ProgressMeter.jl, Interpolations.jl, DifferentialEquations.jl, Suppressor.jl, Printf.jl
```

## Examples
> **Note**
> The document has not been created yet, so please check the examples instead.

`./examples/Fig1_MSLmodel_two_coupled.ipynb`

<img src="https://raw.githubusercontent.com/takyamamoto/GeneralizedPhaseModel.jl/master/figures/fig1.png"> 


## Reference
- Kurebayashi, W., Yamamoto, T., Shirasaka, S., & Nakao, H. (2022). Phase reduction of strongly coupled limit-cycle oscillators. Physical Review Research, 4(4), 043176. <https://doi.org/10.1103/PhysRevResearch.4.043176>
- Kurebayashi, W., Shirasaka, S., & Nakao, H. (2013). Phase reduction method for strongly perturbed limit cycle oscillators. Physical Review Letters, 111(21), 214101. <https://doi.org/10.1103/PhysRevLett.111.214101>

