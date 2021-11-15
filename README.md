# GeneralizedPhaseModel.jl

[![Build Status](https://travis-ci.com/takyamamoto/GeneralizedPhaseModel.jl.svg?branch=main)](https://travis-ci.com/takyamamoto/GeneralizedPhaseModel.jl)
[![Coverage](https://codecov.io/gh/takyamamoto/GeneralizedPhaseModel.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/takyamamoto/GeneralizedPhaseModel.jl)
[![Coverage](https://coveralls.io/repos/github/takyamamoto/GeneralizedPhaseModel.jl/badge.svg?branch=main)](https://coveralls.io/github/takyamamoto/GeneralizedPhaseModel.jl?branch=main)

Julia package for generalized phase reduction method.

<img src="https://raw.githubusercontent.com/takyamamoto/GeneralizedPhaseModel.jl/master/figures/fig1.png"> 

## Requirements
Julia above ver 1.6. 

### Installation
Run the following command. 
```julia
pkg > add https://github.com/takyamamoto/GeneralizedPhaseModel.jl.git
```

### Dependent packages 
```
ForwardDiff.jl, LinearAlgebra.jl, Random.jl, Statistics.jl, ProgressMeter.jl, Interpolations.jl, DifferentialEquations.jl, Suppressor.jl, Printf.jl
```

## Reference
- K. Wataru, T. Yamamoto, S. Shirasaka, and H. Nakao. "Phase reduction of strongly coupled limit-cycle oscillators." *in prep*.
- K. Wataru, S. Shirasaka, and H. Nakao. "Phase reduction method for strongly perturbed limit cycle oscillators." Physical review letters 111.21 (2013): 214101.

