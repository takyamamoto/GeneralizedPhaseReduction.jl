# -*- coding: utf-8 -*-
"""
Takuto Yamamoto
"""
function determining_synchronization(N, D, κ, IΘ, ωI, ξθI, initθgpm, dt, eps=1e-4, Nt=20000)
    """
    Args:
    - N  : 振動子の数
    - D  : 振動子の状態数
    - κ : 振動子の時定数 ∈ R^N

    Return:
    - bool, state?
    """
    Θgpm = zeros(Nt, N)
    integrator_gpm = get_ode_integrator(generalized_coupled_updateΘ, initθgpm, dt, (N, κ, IΘ, ωI, ξθI), RK4())
    for tt in 1:Nt
        # update
        step!(integrator_gpm, dt, true)
        Θgpm[tt, :] = copy(integrator_gpm.u)
        #Θgpm[tt, :] = mod.(copy(integrator_gpm.u), 2π)

        for s in 1:tt-1 # s=t-Δt, t-2Δt, ..., 1
            Θsmt = mod.(Θgpm[s, :] - Θgpm[tt, :], 2π)
            δΘ = mod.(Θgpm[tt+1, :] - Θgpm[tt, :], 2π)
            #Θsmt = Θgpm[s, :] - Θgpm[tt, :]
            #δΘ = Θgpm[tt+1, :] - Θgpm[tt, :]

            # 0をまたぐときの補正
            #Θsmt += (abs.(Θsmt) .> π) .* sign.(Θsmt) * (-2π)
            #δΘ += (abs.(δΘ) .> π) .* sign.(δΘ) * (-2π)

            c = dot(Θsmt, δΘ) / norm(δΘ)
            dist = norm(Θsmt - c*δΘ)
            if 0 <= c <= 1 && dist < eps
                println("c=", c, ", dist=", dist)
                return true, Θgpm
            end
        end

        if tt == Nt - 1
            return false, Θgpm
        end
    end
end

function ArnoldTongue(N, ωI, ξθI, minK=0.05, maxK=0.9, epsk=0.05, Nκ = 100, epsκ = 0.05)
    NK = round(Int, (maxK - minK) / epsk) + 1

    # for recording
    K_list = [] # strength 0.1 ~ 0.9
    κ_list = [] # Δω_list = (κ_list .- 1) * ωI(0)
    synchro_list = [] # true or false
    #dtarray_ = 10 .^(range(-1, -3,length=NK))
    dt = 1e-4

    # init
    κ = ones(N)
    count = 1 # 条件数
    κdirection = 0 # 0 is forward, 1 is backward

    initθgpm = zeros(N) #2pi * rand() .* ones(N) #2pi .* rand(N) [0.0, 0.0]
    for k in 1:NK
        K = (k-1) * epsk + minK # K is connectivity strength 0.1, 0.2, ..., 0.9
        g(X) = G(X, K)
        IΘ = ComputeIsΘ(g, N, D, XsI)
        for j in 1:Nκ
            println("count:", count, ", K=", K, ", κ=", κ[2])
            synchro, _ = determining_synchronization(N, D, κ, IΘ, ωI, ξθI, initθgpm, dt)
            #synchro, θstate = determining_synchronization(N, D, κ, IΘ, ωI, ξθI, initθgpm, dtarray_[k])
            println(synchro) #, θstate)
            # record
            append!(K_list, K) # strength 0.1 ~ 0.9
            append!(κ_list, κ[2]) # Δω_list = (κ_list .- 1) * ωI(0)
            append!(synchro_list, synchro) # true or false

            # update
            count += 1
            if j == 1
                κdirection = synchro ? 0 : 1 # if synchro, κdirection = 0
            end

            if synchro
                if κdirection == 0 κ[2] += epsκ else break end
            else
                if κdirection == 1 κ[2] -= epsκ else break end
            end
        end
    end
    return K_list, κ_list, synchro_list
end

function updateX!(X, I, F, dt, mode="Euler")
    """
    Args:
    - X  : State vector ∈ R^D
    - I  : External force strength; (Scalar)
    - F(X, I)  : Vector field function of the oscillator model; (list of functions)
    - dt : Time step for Euler integration

    Examples:
    > F(X, I) = [dv(X, I), dm(X), dh(X), dn(X)]
    """
    if mode == "RK4"
        k1 = F(X, I) * dt
        k2 = F(X + 0.5*k1, I) * dt
        k3 = F(X + 0.5*k2, I) * dt
        k4 = F(X + k3, I) * dt
        X[:] += (k1 + 2k2 + 2k3 + k4) / 6
    else
        X[:] += F(X, I) * dt
    end
end

function updateXN!(X, I, F, N, dt, mode="Euler")
    """
    Args:
    - X  : State vector ∈ R^D
    - I  : External force strength; (Scalar)
    - F(X, I)  : Vector field function of the oscillator model; (list of functions)
    - N  : number of units
    - dt : Time step for Euler integration

    Examples:
    > F(X, I) = [dv(X, I), dm(X), dh(X), dn(X)]
    """
    for i in 1:N
        if mode == "RK4"
            k1 = F(X[i, :], I[i]) * dt
            k2 = F(X[i, :] + 0.5*k1, I[i]) * dt
            k3 = F(X[i, :] + 0.5*k2, I[i]) * dt
            k4 = F(X[i, :] + k3, I[i]) * dt
            X[i, :] += (k1 + 2k2 + 2k3 + k4) / 6
        else
            X[i, :] += F(X[i, :], I[i]) * dt
        end
    end
end

function Approxω(I, ωI, Imin, Imax, dI)
    """
    Function to calculate the approximate value of ω(I)
    with lookup table style and linear interpolation using ωI.

    Args:
    - θ    : Phase (Float)
    - ω(I) : Angular frequency ∈ R^NI
    - Nθ   : The number of splitting phase for one period; (Int)

    Returns:
    - Approximate ω(I)

    Linear interpolation: f(x) = f(x0) + (x - x0)/(x1 - x0) * (f(x1) - f(x0))
    - ω(I) = ω(I0) + (I - I0)/dI * (ω(x1) - ω(x0))
    - I0 = Imin + dI*(idxf - 1)
    - (I - I0)/dI = (I - Imin - dI*(idxf - 1))/dI = idx - idxf
    """

    @assert (Imin ≤ I ≤ Imax) "Inputted I="*string(I)*" is not in range [Imin, Imax]=["*string(Imin)*", "*string(Imax)*"]"

    # linear interpolation
    idx = (I - Imin)/dI + 1
    idxf = floor(Int, idx)
    idxc = ceil(Int, idx)
    return ωI[idxf] + (idx - idxf)*(ωI[idxc] - ωI[idxf])
end

function Approxζ(θ, I, ζθI, Nθ, Imin, Imax, dI)
    """
    Function to calculate the approximate value of ζ(θ, I)
    with lookup table style using ξθI.

    Args:
    - θ    : Phase (Float)
    - I    : External force strength; (Scalar)
    - ζ(θ, I) : Phase sensitivity function ∈ R^(Nθ x NI)
    - Nθ   : The number of splitting phase for one period; (Int)
    - Imin : Minimum external force strength; (Scalar)
    - Imax : Maximum external force strength; (Scalar)
    - dI   : step of I

    Returns:
    - Approximate ζ(θ, I)
    """

    θ = mod(θ, 2π) # θ ∈[0, 2π]
    θ = θ==2π ? 0 : θ # if θ = 2π, set zero.
    idxθ = floor(Int, Nθ * θ / 2π + 1) # index of θ

    @assert (Imin ≤ I ≤ Imax) "Inputted I="*string(I)*" is not in range [Imin, Imax]=["*string(Imin)*", "*string(Imax)*"]"
    idxI = (I - Imin)/dI + 1
    idxIf = floor(Int, idxI)
    idxIc = ceil(Int, idxI)
    return ζθI[idxθ, idxIf] + (idxI - idxIf)*(ζθI[idxθ, idxIc] - ζθI[idxθ, idxIf])
end

function Approxξ(θ, I, ξθI, Nθ, Imin, Imax, dI)
    """
    Function to calculate the approximate value of ξ(θ, I)
    with lookup table style using ξθI.

    Args:
    - θ    : Phase (Float)
    - I    : External force strength; (Scalar)
    - ξ(θ, I) : Phase sensitivity function ∈ R^(Nθ x NI)
    - Nθ   : The number of splitting phase for one period; (Int)
    - Imin : Minimum external force strength; (Scalar)
    - Imax : Maximum external force strength; (Scalar)
    - dI   : step of I

    Returns:
    - Approximate ξ(θ, I)
    """

    θ = mod(θ, 2π) # θ ∈[0, 2π]
    θ = θ==2π ? 0 : θ # if θ = 2π, set zero.
    idxθ = floor(Int, Nθ * θ / 2π + 1) # index of θ

    @assert (Imin ≤ I ≤ Imax) "Inputted I="*string(I)*" is not in range [Imin, Imax]=["*string(Imin)*", "*string(Imax)*"]"
    idxI = (I - Imin)/dI + 1
    idxIf = floor(Int, idxI)
    idxIc = ceil(Int, idxI)
    return ξθI[idxθ, idxIf] + (idxI - idxIf)*(ξθI[idxθ, idxIc] - ξθI[idxθ, idxIf])
end

function ApproxXs(θ, I, XsI, Nθ, Imin, Imax, dI)
    """
    Function to calculate the approximate value of Xs(θ, I)
    with lookup table style using XsI.

    Args:
    - θ    : Phase (Float)
    - I    : External force strength; (Scalar)
    - Xs(θ, I): stable periodic solution ∈ R^(D x Nθ x NI)
    - Nθ   : The number of splitting phase for one period; (Int)
    - Imin : Minimum external force strength; (Scalar)
    - Imax : Maximum external force strength; (Scalar)
    - dI   : step of I

    Returns:
    - Approximate Xs(θ, I)
    """

    θ = mod(θ, 2π) # θ ∈[0, 2π]
    θ = θ==2π ? 0 : θ # if θ = 2π, set zero.
    idxθ = floor(Int, Nθ * θ / 2π + 1) # index of θ

    @assert (Imin ≤ I ≤ Imax) "Inputted I="*string(I)*" is not in range [Imin, Imax]=["*string(Imin)*", "*string(Imax)*"]"
    idxI = (I - Imin)/dI + 1
    idxIf = floor(Int, idxI)
    idxIc = ceil(Int, idxI)
    return XsI[:, idxθ, idxIf] + (idxI - idxIf)*(XsI[:, idxθ, idxIc] - XsI[:, idxθ, idxIf])
end

function Approxθ(X, I, XsI, Nθ, Imin, Imax, dI)
    """
    Function to calculate the approximate value of θ
    with lookup table style using XsI.

    Args:
    - X    : State variables
    - I    : External force strength; (Scalar)
    - Xs(θ, I): stable periodic solution ∈ R^(D x Nθ x NI)
    - Nθ   : The number of splitting phase for one period; (Int)
    - Imin : Minimum external force strength; (Scalar)
    - Imax : Maximum external force strength; (Scalar)
    - dI   : step of I

    Returns:
    - Approximate θ
    """

    @assert (Imin ≤ I ≤ Imax) "Inputted I="*string(I)*" is not in range [Imin, Imax]=["*string(Imin)*", "*string(Imax)*"]"
    idxI = round(Int, (I - Imin)/dI + 1)
    idxθ = findmin(sum(abs.(XsI[:, :, idxI] .- X), dims=1))[2][2]
    θ = 2pi / Nθ * (idxθ - 1)
    return θ
end

function generalized_coupled_updateΘ!(Θ, IΘ, N, ωI, ξθI, dt)
    XiΘ = zeros(N, N)
    ωΘ = zeros(N)
    Θlim = mod.(Θ, 2π) # θ ∈[0, 2π]

    for i in 1:N
        Ii = IΘ[i](Θlim...)
        ωΘ[i] = ωI(Ii)
        XiΘ[i, :] = - ξθI(Θlim[i], Ii) * Interpolations.gradient(IΘ[i], Θlim...)
    end

    XiΘ += LinearAlgebra.I

    dΘ = inv(XiΘ) * ωΘ # generalized phase equation
    Θ[:] += dΘ * dt
end

flatten(X) = vcat(X...)

function determining_synchronization(N, D, κ, IΘ, ωI, ξθI, initθgpm, dt, eps=1e-3, Nt=20000)
    """
    Args:
    - N  : 振動子の数
    - D  : 振動子の状態数
    - κ : 振動子の時定数 ∈ R^N

    Return:
    - bool, state?
    """
    Xgpm = zeros(Nt, Int(N*D))
    integrator_gpm = get_ode_integrator(generalized_coupled_updateΘ, initθgpm, dt, (N, κ, IΘ, ωI, ξθI), Tsit5())
    for tt in 1:Nt
        # update
        step!(integrator_gpm, dt, true)
        θgpm = mod.(copy(integrator_gpm.u), 2π)
        Xgpm[tt, :] = flatten([XsI[j](θgpm[i], IΘ[i](θgpm...)) for j in 1:D for i in 1:N])

        for s in 1:tt-1 # s=t-Δt, t-2Δt, ..., 1
            Xsmt = Xgpm[s, :] - Xgpm[tt, :]
            δX = Xgpm[tt+1, :] - Xgpm[tt, :]
            c = dot(Xsmt, δX) / norm(δX)
            dist = norm(Xsmt - c*δX)
            if 0 <= c <= 1 && dist < eps
                println("c=", c, ", dist=", dist)
                return true, θgpm
            end
        end

        if tt == Nt - 1
            return false, θgpm
        end
    end
end

function ComputeIsΘ(G, N, D, XsI, NΘ::Int=50, ttmax::Int=50, ϵ=1e-5)
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

    for m in 1:NΘ
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

function generalized_coupled_updateΘ(Θ, p, t)
    N, κ, IΘ, ωI, ξθI = p
    XiΘ = zeros(N, N)
    ωΘ = zeros(N)
    Θlim = mod.(Θ, 2π) # θ ∈[0, 2π]

    for i in 1:N
        Ii = IΘ[i](Θlim...)
        ωΘ[i] = κ[i] * ωI(Ii)
        XiΘ[i, :] = - ξθI(Θlim[i], Ii) * Interpolations.gradient(IΘ[i], Θlim...)
    end

    XiΘ += LinearAlgebra.I
    dΘ = XiΘ \ ωΘ # generalized phase equation
    return dΘ
end

state(θ, Xs) = Xs[:, floor(Int, mod(θ, 2π)/2π * size(Xs)[2])+1];

function ArnoldTongueOriginalSystem(
        func, connect_func, dims, initX0, dt, alg;
        eps=1e-8, λthr=1e-6, exploreDir=0, initK=1e-2, endK=0.5, Δk=1e-2, Δκ = 1e-4, maxNκ=200,
        print_progress=false)
    """
    Args
    - connect_func : connectivity function
    - initX0：
    - eps   : convergence check threshod
    - λthr : eigen value threshold
    - exploreDir：exploration direction: 0 is forward, 1 is backward

    Returns
    - listK :strength 0.1 ~ 0.9
    - listκ : Δω_list = (listκ .- 1) * ωI(0)
    - listStabilities : boolean
    - listXt :
    """

    # initial variables
    listK, listκ, listStabilities, listXt = [], [], [], [] # for recording
    NK = round(Int, (endK - initK) / Δk) + 1
    κ = ones(N) #
    count = 1 # condition number

    @showprogress for k in 1:NK
        K = (k-1) * Δk + initK # K is connectivity strength 0.1, 0.2, ...,
        g(X) = connect_func(X, K)
        for j in 1:maxNκ
            if print_progress println("count:", count, ", K=", K, ", κ=", κ[2]) end;
            input_params = [g, κ]
            if k == 1 && j == 1
                global initT, _, Xs = FindStablePeriodicSolution(func, input_params, dims, Int(1/dt), initX0, dt, alg, print_progress=false);
            end
            convergence, stability, s, T,  Xt = ShootingMethod(func, input_params, dims, Xs, initT, dt, alg, eps, λthr, print_progress=false);
            Xtu = convergence ? cat(Xt.u..., dims=2) : nothing
            if print_progress 
                println("stability: ", stability, ", size(Xs): ", size(Xs)) 
                println("s: ", s)
            end;

            # recording
            append!(listK, K); append!(listκ, κ[2]); append!(listStabilities, stability); append!(listXt, [Xtu])

            # update
            count += 1
            if j == 1 exploreDir = stability ? 0 : 1 end # if stability, exploreDir = 0
            if stability
                global initT, Xs = T, Xtu # update initT and Xs
                if exploreDir == 0 κ[2] += Δκ; else break; end
            else
                if exploreDir == 1 κ[2] -= Δκ; else break; end
            end
            @assert j < maxNκ "Divergence occurred."
        end
    end
    return listK, listκ, listStabilities, listXt
end

function FindArnoldTongueBorder(listK, listκ, listStabilities)
    uniqueK = unique(listK)
    Δω = (listκ .- 1)
    NK = size(uniqueK)[1]
    bordersidx = zeros(NK)
    bordersΔω = zeros(NK)
    cumsumidx = 0
    for l in 1:NK
        k = uniqueK[l]
        idxk = (listK .== k)
        stability_k = listStabilities[idxk]
        Δωk = Δω[idxk]
        #println(l, stability_k, sum(stability_k))
        if size(stability_k)[1] == 2
            idxkborder = 1
            border_k = sum(Δωk)/2
        else
            n_one = sum(stability_k)
            n_zero = size(stability_k)[1] - n_one
            if n_one >= n_zero
                idxkborder = findall(x->x==true, stability_k)[end]
            else
                idxkborder = findall(x->x==false, stability_k)[end]
            end
            border_k = (Δωk[idxkborder] + Δωk[idxkborder+1])/2
        end
        bordersidx[l], bordersΔω[l] = cumsumidx+idxkborder, border_k
        cumsumidx += size(stability_k)[1]
    end
    return Int.(bordersidx), bordersΔω
end

function shooting_method(
        func, input_param, dims, initX0, initT, dt, alg, eps, λthr, maxiter=100;
        print_progress=true)

    """
    func : vector field
    input_param : G (connectivity matrix) or I (input) etc...
    dims : dimensions of system.
    initX0 : initial X0
    initT : initial period T
    eps : small value for convergence
    λthr : eigenvalue thr
    maxiter : maximum iteration
    plot : bool. plot X or not.
    """
    Jfunc(X, input_param) = ForwardDiff.jacobian(X -> func(X, input_param...), X)

    function updateTfunc!(dX, X, p, t)
        input_param, T = p
        dX[:] = T*func(X, input_param...)
    end

    function updateY!(dY, Y, p, t)
        input_param, Xt, T = p
        J = [T*Jfunc(Xt(t), input_param) func(Xt(t), input_param...); zeros(1, size(Y)[1])]
        dY[:, :] = J * Y
    end

    # set initial states
    X0 = initX0
    T = initT
    for iter = 1:maxiter
        Xt = get_ode_solution(updateTfunc!, X0, (0.0, 1.0), dt, (input_param, T), alg)
        ΔX = norm(Xt.u[end, :][1] - X0)
        if print_progress println("[$iter] : X(1)-X(0) = $(@sprintf("%.3E", ΔX)), T= $(@sprintf("%.3f", T))") end;

        # Convergence
        if ΔX < eps
            global convergence, Xtc, X0c, Tc  = true, Xt, X0, T
            break
        end

        Yt = get_ode_solution(updateY!, eye(dims+1), (0.0, 1.0), dt, (input_param, Xt, T), alg, save_everystep=false)

        # Update (Newton method)
        U, s, V = svd(Yt.u[end] - LinearAlgebra.I)

        idx = argmin(abs.(s))
        s2 = 1.0 ./ (s .+ 1e-8) # avoid zero division
        s2[maximum(s)*λthr .>= s] .= 0
        #s2[idx] .= 0
        dX = V * Diagonal(s2) * U' * vcat(Xt[end] - X0, 0)
        X0[:] -= dX[1:end-1]
        T -= dX[end]

        if T <= dt || iter == maxiter
            global convergence, Xtc, X0c, Tc = false, nothing, nothing, nothing
            break
        end
    end
    return convergence, Xtc, X0c, Tc
end

function JudgeStability(func, input_param, dims, T, X0, dt, alg)
    """
    F : vector field
    input_param : G (connectivity matrix) or I (input) etc...
    T : period
    dims : dimensions of state
    """

    Jfunc(X, input_param) = ForwardDiff.jacobian(X -> func(X, input_param...), X)

    function updatefunc!(dX, X, input_param, t)
        dX[:] = func(X, input_param...)
    end

    function updateYs!(dY, Y, p, t)
        input_param, Xt = p
        dY[:] = Jfunc(Xt(t), input_param) * Y
    end

    Xt = get_ode_solution(updatefunc!, X0, (0.0, T), dt, input_param, alg)
    Yt = get_ode_solution(updateYs!, eye(dims), (0.0, T), dt, (input_param, Xt), alg, save_everystep=false)
    _, s, _ = svd(Yt.u[end])
    idx = argmin(abs.(s .- 1.0))
    return all(s[1:end .!= idx] .< 1.0), s
end

function ShootingMethod(
        func, input_param, dims, Xs, initT, dt, alg, eps, λthr; max_trial=5, print_progress=false)
    """
    dt  : e.g. 1e-3
    alg : e.g. Tsit5()
    e.g.
    # Generate initial states and period.
    initX = randn(dims)
    initT, _, Xs = FindStablePeriodicSolution(func, input_param, dims, Int(1/dt), initX, dt, alg);
    """

    # Try shooting method until convergence (maximum number of times is `maxtrial`)
    initX0 = Xs[:, end];
    @suppress_err global convergence, Xt, X0, T = try
        shooting_method(func, input_param, dims, initX0, initT, dt, alg, eps, λthr, print_progress=print_progress)
    catch
        false, nothing, nothing, nothing
    end

    # Check stability
    if convergence
        stability, s = JudgeStability(func, input_param, dims, T, X0, dt, alg)
    else
        stability = false
        s = nothing
    end
    return convergence, stability, s, T, Xt
end
