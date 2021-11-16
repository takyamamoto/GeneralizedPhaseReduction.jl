module GeneralizedPhaseReduction
    export eye, vec, 
           update_z!, find_stable_periodic_solution, phase_sensitivity_func, approx_z,
           generalized_phase_sensitivity_func,
           get_ode_integrator, get_ode_solution,
           compute_QΘ, compute_IsΘ, conventinal_coupled_updateΘ, 
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
    using Printf

    include("utils.jl")
    include("phase_sensitivity.jl")
    include("generalized_phase_sensitivity.jl")
    include("differential_equations_extensions.jl")
    include("coupled_phase_equation.jl")
end
