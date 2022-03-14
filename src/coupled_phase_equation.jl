# -*- coding: utf-8 -*-
function compute_QΘ(G, N, D, XsI, κ, ωI, NΘ::Int=50, ttmax::Int=50, ϵ=1e-5, ρ=0.1, dt=1e-2)
    """
    Compute Q(θ)

    Args:
    - N = 2 # number of units
    - D   : Dimension of states; (Int)
    - Nθ  : The number of splitting phase for one period; (Int)

    Returns:
    - Q(θ)   :  ∈ R^N x NΘ x NΘ
    """

    X = zeros(N, D)
    Θ, Θ̃ = zeros(N), zeros(N)
    Θrange = range(0, 2π, length=NΘ)
    Qₜ = zeros(N, NΘ, NΘ)
    α = ρ*dt
    
    @showprogress "Computing P(θ₁, θ₂)..." for tt in 1:ttmax
        Qₜ₋₁ = copy(Qₜ)
        
        for m in 1:NΘ
            for n in 1:NΘ
                Θ = [Θrange[m], Θrange[n]]
                for s in 1:-dt:0 
                    for i in 1:N
                        Θ̃[i] = mod(Θ[i] - κ[i] * ωI(Qₜ[i, m, n]) * s, 2π)
                    end
                    #Θ̃indices = clamp.(floor.(Int, Θ̃/(2π)*NΘ) .+ 1, 1, NΘ)
                    for i in 1:N
                        X[i, :] = [XsI[j](Θ̃[i], Qₜ[i, m, n]) for j in 1:D] # for approximation
                        #X[i, :] = [XsI[j](Θ̃[i], Qₜ[i, Θ̃indices...]) for j in 1:D]
                    end
                    Qₜ[:, m, n] += α*(G(X) - Qₜ[:, m, n]) # low-pass filter
                end
            end
        end
        if sum(abs.(Qₜ .- Qₜ₋₁))/(NΘ^2) < ϵ # check convergence
            break
        end
    end
    sitpQΘ = [scale(interpolate(Qₜ[i, :, :], BSpline(Cubic(Line(OnGrid())))), Θrange, Θrange) for i in 1:N]
    return sitpQΘ
end

function compute_IsΘ(G, N, D, XsI, NΘ::Int=50, ttmax::Int=50, ϵ=1e-5)
    """
    Compute I₀(θ)

    Args:
    - N = 2 # number of units
    - D   : Dimension of states; (Int)
    - Nθ  : The number of splitting phase for one period; (Int)

    Returns:
    - I₀(θ)   :  ∈ R^N x NΘ x NΘ
    """

    X = zeros(N, D)
    Θ = zeros(N)
    Θrange = range(0, 2π, length=NΘ)

    IΘ = zeros(N, NΘ, NΘ)

    @showprogress "Computing I₀(θ₁, θ₂)..." for m in 1:NΘ
        for n in 1:NΘ
            Θ = [Θrange[m], Θrange[n]]
            # initial I
            if n > 1
                Iₜ = copy(IΘ[:, m, n-1])
            else
                Iₜ = zeros(N)
            end

            for tt in 1:ttmax
                for i in 1:N
                    X[i, :] = [XsI[j](Θ[i], Iₜ[i]) for j in 1:D]
                end

                Iₜ₋₁ = copy(Iₜ)
                Iₜ = G(X) # update

                if sum(abs.(Iₜ .- Iₜ₋₁)) <  ϵ # check convergence
                    break
                end
            end

            IΘ[:, m, n] = Iₜ # save
        end
    end

    sitpIΘ = [scale(interpolate(IΘ[i, :, :], BSpline(Cubic(Line(OnGrid())))), Θrange, Θrange) for i in 1:N]
    return sitpIΘ
end

function conventinal_coupled_updateΘ(Θ, p, t)
    N, D, κ, G, ωI, ζθI, XsI = p
    return κ * ωI(0) + ζθI.(mod.(Θ, 2π), 0) .* G(hcat([[XsI[j](mod(θ, 2π), 0) for j in 1:D] for θ in Θ]...)')
end

function generalized_coupled_updateΘ_I(Θ, p, t)
    N, κ, IΘ, ωI, ξθI = p
    Θ = mod.(Θ, 2π) # θ ∈[0, 2π]
    inp = [IΘ[i](Θ...) for i in 1:N]
    ωΘ = [κ[i] * ωI(inp[i]) for i in 1:N]
    XiΘ = LinearAlgebra.I　- hcat([ξθI(Θ[i], inp[i]) * Interpolations.gradient(IΘ[i], Θ...) for i in 1:N]...)'
    dΘ = XiΘ \ ωΘ # generalized phase equation
    return dΘ
end

function generalized_coupled_updateΘ_PQ(Θ, p, t)
    N, D, κ, G, QΘ, ωI, ζθI, ξθI, XsI = p
    Θ = mod.(Θ, 2π) # θ ∈[0, 2π]
    q̂ = [QΘ[i](Θ...) for i in 1:N]
    p̂ = G(hcat([[XsI[j](Θ[i], q̂[i]) for j in 1:D] for i in 1:N]...)') - q̂
    ϕ = [κ[i] * ωI(q̂[i]) + ζθI(Θ[i], q̂[i])*p̂[i] for i in 1:N]
    XiΘ = LinearAlgebra.I - hcat([ξθI(Θ[i], q̂[i]) * Interpolations.gradient(QΘ[i], Θ...) for i in 1:N]...)'
    dΘ = XiΘ \ ϕ # generalized phase equation
    return dΘ
end

function phase2cum_phase(Θ, N)
    ΔΘ = Θ[2:end, :] - Θ[1:end-1, :]
    idx = ΔΘ .< -1.9π # detect phase jump point
    ΔΘ[idx] .= mod.(ΔΘ[idx], 2π);
    Θcumsum = cumsum(ΔΘ, dims=1);
    for i in 1:N
        Θcumsum[:, i] .+= Θ[1, i]
    end
    return Θcumsum = [Θ[1, :]'; Θcumsum]
end

# +
# for fixed time steps
function coupled_original_system(N::Int, D::Int, Nt::Int, dt::Float64, XsI, G, coupled_func, initθ, κ, alg=Tsit5())
    X = zeros(Nt, N, D)
    Θg, Θc = zeros(Nt, N), zeros(Nt, N);
    initX = hcat([[XsI[j](mod(θ, 2π), 0) for j in 1:D] for θ in initθ]...)'
    integrator = get_ode_integrator(coupled_func, initX, dt, (G, κ), alg)
    for tt in 1:Nt
        x = copy(integrator.u)
        X[tt, :, :] = x # memory
        Θg[tt, :] = mod.(atan.(x[:, 2], x[:, 1] - G(x)), 2π)
        Θc[tt, :] = mod.(atan.(x[:, 2], x[:, 1]), 2π)
        step!(integrator, dt, true)
    end
    Θg_cumsum = phase2cum_phase(Θg, N)
    Θc_cumsum = phase2cum_phase(Θc, N)
    return X, Θg_cumsum, Θc_cumsum
end

# for adaptive time step
function coupled_original_system(N::Int, D::Int, T::Float64, XsI, G, coupled_func, initθ, κ, alg=Tsit5())
    initX = hcat([[XsI[j](mod(θ, 2π), 0) for j in 1:D] for θ in initθ]...)'
    sol = get_ode_solution(coupled_func, initX, (0, T), nothing, (G, κ), alg)
    Nt = length(sol.u)
    X = zeros(Nt, N, D)
    Θg, Θc = zeros(Nt, N), zeros(Nt, N);
    for tt in 1:Nt
        x = sol.u[tt]
        X[tt, :, :] = x # memory
        Θg[tt, :] = mod.(atan.(x[:, 2], x[:, 1] - G(x)), 2π)
        Θc[tt, :] = mod.(atan.(x[:, 2], x[:, 1]), 2π)
    end
    Θg_cumsum = phase2cum_phase(Θg, N)
    Θc_cumsum = phase2cum_phase(Θc, N)
    
    return sol.t, X, Θg_cumsum, Θc_cumsum
end

# +
# for fixed time steps
function coupled_conventinal_phase_model(N::Int, D::Int, Nt::Int, dt::Float64, XsI, G, ωI, ζθI, initθ, κ, alg=Tsit5())
    X = zeros(Nt, N, D) # states
    Θ = zeros(Nt, N)    # phase
    integrator = get_ode_integrator(conventinal_coupled_updateΘ, initθ, dt, (N, D, κ, G, ωI, ζθI, XsI), alg)
    for tt in 1:Nt
        θ = mod.(copy(integrator.u), 2π)
        Θ[tt, :] = copy(integrator.u)
        for i in 1:N
            X[tt, i, :] = [XsI[j](θ[i], 0) for j in 1:D]
        end
        step!(integrator, dt, true) # update
    end
    return X, Θ
end

# for adaptive time step
function coupled_conventinal_phase_model(N::Int, D::Int, T::Float64, XsI, G, ωI, ζθI, initθ, κ, alg=Tsit5())
    sol = get_ode_solution(conventinal_coupled_updateΘ, initθ, (0, T), nothing, (N, D, κ, G, ωI, ζθI, XsI), alg);
    Θ = hcat(sol.u...)';
    Nt = length(sol.u)
    X = zeros(Nt, N, D)
    for tt in 1:Nt
        for i in 1:N
            X[tt, i, :] = [XsI[j](mod(Θ[tt, i], 2π), 0) for j in 1:D] # memory
        end
    end    
    return sol.t, X, Θ
end

# +
# for fixed time steps
function coupled_generalized_phase_model_I(N::Int, D::Int, Nt::Int, dt::Float64, XsI, IΘ, ωI, ξθI, initθ, κ, alg=Tsit5())
    X = zeros(Nt, N, D) # states
    Θ = zeros(Nt, N)    # phase
    integrator = get_ode_integrator(generalized_coupled_updateΘ_I, initθ, dt, (N, κ, IΘ, ωI, ξθI), alg)
    for tt in 1:Nt
        θ = mod.(copy(integrator.u), 2π)
        Θ[tt, :] = copy(integrator.u)
        for i in 1:N
            X[tt, i, :] = [XsI[j](θ[i], IΘ[i](θ...)) for j in 1:D]
        end
        step!(integrator, dt, true) # update
    end
    return X, Θ
end

# for adaptive time step
function coupled_generalized_phase_model_I(N::Int, D::Int, T::Float64, XsI, IΘ, ωI, ξθI, initθ, κ, alg=Tsit5())
    sol = get_ode_solution(generalized_coupled_updateΘ_I, initθ, (0, T), nothing, (N, κ, IΘ, ωI, ξθI), alg);
    Θ = hcat(sol.u...)';
    Nt = length(sol.u)
    X = zeros(Nt, N, D)
    for tt in 1:Nt
        for i in 1:N
            X[tt, i, :] = [XsI[j](mod(Θ[tt, i], 2π), IΘ[i](mod.(Θ[tt, :], 2π)...)) for j in 1:D] # memory
        end
    end    
    return sol.t, X, Θ
end

# +
# for fixed time steps
function coupled_generalized_phase_model_PQ(N::Int, D::Int, Nt::Int, dt::Float64, XsI, QΘ, ωI, ζθI, ξθI, initθ, κ, G, alg=Tsit5())
    X = zeros(Nt, N, D) # states
    Θ = zeros(Nt, N)    # phase
    integrator = get_ode_integrator(generalized_coupled_updateΘ_PQ, initθ, dt, (N, D, κ, G, QΘ, ωI, ζθI, ξθI, XsI), alg)
    for tt in 1:Nt
        θ = mod.(copy(integrator.u), 2π)
        Θ[tt, :] = copy(integrator.u)
        q̂ = [QΘ[i](θ...) for i in 1:N]
        Iext = G(hcat([[XsI[j](θ[i], q̂[i]) for j in 1:D] for i in 1:N]...)')
        for i in 1:N
            X[tt, i, :] = [XsI[j](θ[i], Iext[i]) for j in 1:D]
        end
        step!(integrator, dt, true) # update
    end
    return X, Θ
end

# for adaptive time step
function coupled_generalized_phase_model_PQ(N::Int, D::Int, T::Float64, XsI, QΘ, ωI, ζθI, ξθI, initθ, κ, G, alg=Tsit5())
    sol = get_ode_solution(generalized_coupled_updateΘ_PQ, initθ, (0, T), nothing, (N, D, κ, G, QΘ, ωI, ζθI, ξθI, XsI), 
                           alg);
    Θ = hcat(sol.u...)';
    Nt = length(sol.u)
    X = zeros(Nt, N, D)
    for tt in 1:Nt
        θ = mod.(Θ[tt, :], 2π)
        q̂ = [QΘ[i](θ...) for i in 1:N]
        Iext = G(hcat([[XsI[j](θ[i], q̂[i]) for j in 1:D] for i in 1:N]...)')
        for i in 1:N
            X[tt, i, :] = [XsI[j](θ[i], Iext[i]) for j in 1:D]
        end
    end    
    return sol.t, X, Θ
end
