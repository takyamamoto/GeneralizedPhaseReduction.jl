module GeneralizedPhaseModel
    #export GPM
    #const GPM = GeneralizedPhaseModel
    export eye, vec, 
           updateZ!, FindStablePeriodicSolution, ComputePhaseSensitivityFunction, ApproxZ,
           ComputeGeneralizedPhaseSensitivityFunction,
           get_ode_integrator, get_ode_solution,
           ComputeQΘ, ComputeIsΘ, low_pass_filter, conventinal_coupled_updateΘ, 
           generalized_coupled_updateΘ_I, generalized_coupled_updateΘ_PQ,
           phase2cum_phase, coupled_original_system, coupled_conventinal_phase_model, 
           coupled_generalized_phase_model_I, coupled_generalized_phase_model_PQ

    using ForwardDiff
    using LinearAlgebra
    using Random
    using Statistics
    using ProgressMeter
    using Interpolations
    using DifferentialEquations
    using Suppressor
    using PyPlot
    using Printf

    include("utils.jl")
    include("phase_sensitivity.jl")
    include("generalized_phase_sensitivity.jl")
    include("differential_equations_extensions.jl")
    include("coupled_phase_equation.jl")
end
