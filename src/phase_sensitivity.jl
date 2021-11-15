# -*- coding: utf-8 -*-
# Adjoint equation with Euler integration
function update_z!(X, Z, input_param, F, JF, ω, dt)
    """
    Args:
    - X  : State vector ∈ R^D
    - Z  : Phase sensitivity vector ∈ R^D
    - input_param  : External force strength; (Scalar)
    - F(X, I)  : Vector field function of the oscillator model; (list of functions)
    - JF(X, I) : Jacobian of F (see the following Examples).
    - ω  : Angular frequency of stable periodic solutions (Hz)
    - dt : Time step for Euler integration

    Examples:
    > F(X, input_param) = [dv(X, I), dm(X), dh(X), dn(X)]
    > JF(X, input_param) = ForwardDiff.jacobian(X -> F(X, I), X)
    """
    # Normalize Z
    s = F(X, input_param...)' * Z
    Z[:] *= ω/s

    # Evolve Z backward
    Z[:] += JF(X, input_param...)' * Z * dt
end

function find_stable_periodic_solution(
    F, input_param, D::Int, Nθ::Int, initX=nothing, dt=1e-4, alg=Tsit5(),
    origin_val_idx::Int=1, origin_thr=nothing;
    TREL=200.0, RMAX=Int(5e6), NTRIAL=3, print_progress=true)

    """
    Find a stable periodic solution and compute period T and frequency ω.

    Args:
    - F(X, input_param) : Vector field function of the oscillator model; (list of functions)
    - I   : External force strength; (Scalar)
    - D   : Dimension of states; (Int)
    - Nθ  : The number of splitting phase for one period; (Int)
    - initX   : initial X states
    - origin_val_idx: index of variable which is criterion for zero phase (default:1) 1≦origin_val_idx≦D
    - origin_thr : threshold of zero phase (default:nothing->mean of X1)
    - dt : Time step for numerical integration
    - nimode : Numerical integration mode ("Euler", "RK4")
    - TREL : Initial relaxation time
    - TWAIT : Waiting time for measuring the oscillation period
    - RMAX : Buffer length

    Returns:
    - T : Period of stable periodic solutions (sec)
    - ω : Angular frequency of stable periodic solutions (Hz)
    - Xs(θ) : Stable periodic solutions ∈ R^(D x N_θ)
    """

    # Define the initial state of the state vector.
    if initX == nothing
        initX = rand(D)
    end

    """Determining the points that correspond to phase zero."""

    @assert (1 <= origin_val_idx <= D) "origin_val_idx must be in range [1, D]."

    function update!(dX, X, input_param, t)
        dX[:] = F(X, input_param...)
    end

    # Discard initial relaxation to the limit-cycle orbit
    if print_progress println("[1/4] Relaxation...") end
    integrator = get_ode_integrator(update!, initX, dt, input_param, alg)

    step!(integrator, TREL, true)

    # If origin_thr is not defined, the state after TREL (sec) from the initial state is set to phase zero.
    if origin_thr == nothing
        NTREL = Int(TREL/dt)
        Xrec = zeros(D, NTREL) # recording array
        for tt in 1:NTREL
            Xrec[:, tt] = integrator.u
            step!(integrator, dt, true)
        end

        origin_thr = median(Xrec[origin_val_idx, :])
    end

    if print_progress println("[2/4] Looking for the origin of the state X...") end
    for trial in 1:NTRIAL
        if trial > 1
            integrator = get_ode_integrator(update!, initX, dt, input_param, alg)
        end
        # Xarr = zeros(RMAX, D)
        for tt in 1:RMAX
            step!(integrator, dt, true)
            Xtm1 = integrator.uprev
            Xt = integrator.u
            #Xarr[tt, :] = Xt
            if tt > 1 && (Xtm1[origin_val_idx] < origin_thr) && (Xt[origin_val_idx] > origin_thr)
                # Linear interpolation: f(x) = f(x0) + (x - x0)/(x1 - x0) * (f(x1) - f(x0))
                global originX = Xtm1 + (origin_thr - Xtm1[origin_val_idx])/(Xt[origin_val_idx] - Xtm1[origin_val_idx]) * (Xt - Xtm1)
                @goto label #break
            elseif tt >= RMAX
                if trial >= NTRIAL
                    break
                end
                # @assert false "Insufficient buffer or wrong time step."
                println("Insufficient time step. Set dt=", dt, "->", dt*10)
                dt *= 10.0 # Multiply the time step by ten. Set the first time step small.)
            end
        end
        if trial >= NTRIAL
            @assert false "Insufficient buffer or wrong time step."
        end
    end

    @label label

    """Numerically integrate with setting the phase of state originX as zero
    and compute until the phase comes to zero (θ(t)=ω(t+t0)=2πk (k ∈ N)) again.
    Calculate the period T, the angular frequency ω = 2π/T．"""

    if print_progress println("[3/4] Measuring the oscillation period...") end

    T = 0 # period
    ω = 0 # angular frequency

    reinit!(integrator, originX)
    #  Measure the oscillation period. Here, X is originX.
    for tt in 1:RMAX
        step!(integrator, dt, true)
        Xtm1 = integrator.uprev
        Xt = integrator.u

        if tt > 1 && (Xtm1[origin_val_idx] < origin_thr) && (Xt[origin_val_idx] > origin_thr)
            T = tt*dt
            ω = 2*pi/T
            break
        elseif (tt >= RMAX)
            @assert false "Insufficient buffer or wrong time step."
        end
    end

    """Perform numerical integration for one period with time step size dt=T/Nθ
    and compute the solution as a stable periodic solution."""

    if print_progress println("[4/4] Computing the stable periodic solution...") end
    Xs = zeros(D, Nθ)
    dtθ = T / Nθ

    # update integrator
    integrator = get_ode_integrator(update!, initX, dtθ, input_param, alg)
    reinit!(integrator, originX)

    # Record X for one period
    for tt in 1:Int(Nθ)
        Xs[:, tt] = integrator.u
        step!(integrator, dtθ, true)
    end

    return T, ω, Xs
end

function phase_sensitivity_func(F, I, D::Int, Nθ::Int, T, Xs::Array, itp="cubic", REP::Int=20)
    """
    Compute phase sensitivity function Z(θ) with the adjoint method.

    Args:
    - F(X, I) : Vector field function of the oscillator model; (list of functions)
    - I   : External force strength; (Scalar)
    - D   : Dimension of states; (Int)
    - Nθ  : The number of splitting phase for one period; (Int)
    - T   : Period of stable periodic solutions (sec); (Float)
    - Xs  : A stable solution of the oscillator model; (Array ∈ R^(D x Nθ))
    - itp : interpolation mode. "linear" or "cubic"
    - REP : Relaxation time for the adjoint equation; (int)

    Returns:
    - Z(θ, I::const.) : Phase sensitivity function ∈ R^(D x Nθ)
    """

    dtθ = T / Nθ
    ω = 2*pi/T
    JF(X, I) = ForwardDiff.jacobian(X -> F(X, I), X)

    Zθ = zeros(D, Nθ) # phase sensitivity function

    Z = ones(D) # Initial phase sensitivity

    for rep in 1:REP
        for tt in Nθ:-1:1
            X = Xs[: ,tt]
            update_z!(X, Z, I, F, JF, ω, dtθ)

            if rep == REP
                Zθ[:, tt] = Z
            end
        end
    end

    """ # if use interpolation, we have to change ComputeGeneralizedPhaseSensitivityFunction.
    if itp == "cubic"
        sitpZθ = [scale(interpolate(Zθ[i, :], BSpline(Cubic(Line(OnGrid())))), θrange) for i in 1:D]
    else
        sitpZθ = [scale(interpolate(Zθ[i, :], BSpline(Linear())), θrange) for i in 1:D]
    end
    """

    return Zθ
end

function approx_z(θ, Zθ, Nθ)
    """
    Function to calculate the approximate value of Z(θ)
    with lookup table style using Zθ.

    Args:
    - θ : Phase (Float)
    - Z(θ, I::const.) : Phase sensitivity function ∈ R^(D x Nθ)
    - Nθ  : The number of splitting phase for one period; (Int)
            Nθ = size(Zθ)[2]

    Returns:
    - Approximate Z(θ, I)
    """

    θ = mod(θ, 2π) # θ ∈[0, 2π]
    idx = round(Int, Nθ * θ / 2π + 1) # index of θ
    return Zθ[:, idx]
end
