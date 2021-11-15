# -*- coding: utf-8 -*-
function generalized_phase_sensitivity_func(
    F, Imin, Imax, dI, dims::Int, Nθ::Int, initX=nothing, dt=1e-4, alg=Tsit5(), 
    origin_val_idx::Int=1, origin_thr=nothing, itp="cubic", extrap=false,
    print_progress=true)
    """
    Compute angular frequency ω(I)
    and generalized phase sensitivity function ζ(θ, I), ξ(θ, I) with the adjoint method.

    Args:
    - F(X, I) : Vector field function of the oscillator model; (list of functions)
    - Imin : Minimum external force strength; (Scalar)
    - Imax : Maximum external force strength; (Scalar)
    - dI  : step of I
    - dims   : Dimension of states; (Int)
    - Nθ  : The number of splitting phase for one period; (Int)
    - initX : Array of initial ∈ R^(D x NI)
    - itp : interpolation mode. "linear" or "cubic"

    Returns: Interpolated object
    - ω(I)    : Angular frequency ∈ R^NI
    - ζ(θ, I) : Phase sensitivity function ∈ R^(Nθ x NI)
    - ξ(θ, I) : Phase sensitivity function ∈ R^(Nθ x NI)
    - Xs(θ, I): stable periodic solution ∈ R^(D x Nθ x NI)
    """
    Irange = Imin:dI:Imax
    NI = size(Irange)[1]

    flag_arraydt = dt isa Array
    @assert (flag_arraydt || (dt isa Number)) "'dt' must be a scalar or an array."
    if flag_arraydt
        @assert size(dt)[1] == NI "The size of the array of 'dt' must be the same as the number of conditions for the external force 'I'."
        dtarray = copy(dt)
    end

    if initX == nothing
        initX = rand(dims)
    end

    flag_initX = size(initX) == (NI, dims)
    if flag_initX
        initXarray = copy(initX)
    else
        @assert size(initX) == (dims, ) "The size of 'initX' must be the same as (dims, )."
    end

    θrange = range(0, 2π, length=Nθ)

    ZθI = zeros(dims, Nθ, NI)
    ωI = zeros(NI)
    XsI = zeros(dims, Nθ, NI)
    @showprogress "[1/3] Computing Xs(θ, I) and Z(θ, I)..." for i in 1:NI
        if flag_arraydt
            dt = dtarray[i]
        end

        if flag_initX
            initX = initXarray[i, :]
        end

        T, ω, Xs = find_stable_periodic_solution(F, Irange[i], dims, Nθ, initX, dt, alg, origin_val_idx, origin_thr, print_progress=false)
        Zθ = phase_sensitivity_func(F, Irange[i], dims, Nθ, T, Xs)

        ωI[i] = ω
        ZθI[:, :, i] = Zθ
        XsI[:, :, i] = Xs
    end

    ζθI = zeros(Nθ, NI)
    ξθI = zeros(Nθ, NI)
    JFI(X, I) = ForwardDiff.derivative(I -> F(X, I), I)

    # Compute ζ(θ, I)
    if print_progress println("[2/3] Computing ζ(θ, I)...") end

    for i in 1:NI
        for th in 1:Nθ
            ζθI[th, i] = JFI(XsI[:, th, i], Irange[i])' * ZθI[:, th, i]
        end
    end

    # numerical integration
    dθ = 2π / Nθ
    ζI = sum(ζθI, dims=1)[1, :] / Nθ # dθ/2π = 1/Nθ

    # scaled interpolated variables
    if itp == "cubic"
        sitpXsI = [scale(interpolate(XsI[i, :, :], BSpline(Cubic(Line(OnGrid())))), θrange, Irange) for i in 1:dims]
    else
        sitpXsI = [scale(interpolate(XsI[i, :, :], BSpline(Linear())), θrange, Irange) for i in 1:dims]
    end

    # Compute ξ(θ, I)
    if print_progress println("[3/3] Computing ξ(θ, I)...") end
    for i in 1:NI
        ξθI[1, i] = - [Interpolations.gradient(sitpXsI[j], 0, Irange[i])[2] for j in 1:dims]' * ZθI[:, 1, i]
    end

    for i in 1:NI
        for th in 2:Nθ
            ξθI[th, i] = ξθI[th-1, i] - (ζθI[th, i] - ζI[i]) * dθ / ωI[i]
        end
    end

    # scaled interpolated variables
    if itp == "cubic"
        sitpωI = scale(interpolate(ωI, BSpline(Cubic(Line(OnGrid())))), Irange)
        sitpζθI = scale(interpolate(ζθI, BSpline(Cubic(Line(OnGrid())))), θrange, Irange)
        sitpξθI = scale(interpolate(ξθI, BSpline(Cubic(Line(OnGrid())))), θrange, Irange)
    else
        sitpωI = scale(interpolate(ωI, BSpline(Linear())), Irange)
        sitpζθI = scale(interpolate(ζθI, BSpline(Linear())), θrange, Irange)
        sitpξθI = scale(interpolate(ξθI, BSpline(Linear())), θrange, Irange)
    end

    if extrap
        etpωI = extrapolate(sitpωI, Flat()) 
        etpζθI = extrapolate(sitpζθI, Flat()) 
        etpξθI = extrapolate(sitpξθI, Flat()) 
        etpXsI = [extrapolate(sitpXsI[i], Flat()) for i in 1:dims] 
        return etpωI, etpζθI, etpξθI, etpXsI
    else
        return sitpωI, sitpζθI, sitpξθI, sitpXsI
    end
end
