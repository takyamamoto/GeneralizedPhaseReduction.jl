function get_ode_integrator(func, u0, dt=nothing, param=nothing, alg=Tsit5(), maxiters=Inf;
    reltol=1e-12, abstol=1e-12)
    """
    Args:
    - func : the update function
    - u0  : Inital vector of value
    - param   : additional param
    - dt  : Time step for integration
    - alg : update algorithm (e.g. Euler(), RK4())
    - maxiters : maxiters of integration. Default is Inf.
    """

    tspan = (0.0, Inf)

    if param == nothing
        prob = ODEProblem(func, u0, tspan)
    else
        prob = ODEProblem(func, u0, tspan, param)
    end

    if dt == nothing
        integrator = init(prob, alg, adaptive=true, maxiters=maxiters)
    else
        integrator = init(prob, alg, dt=dt, adaptive=false, maxiters=maxiters)
    end
    integrator.opts.reltol = reltol
    integrator.opts.abstol = abstol
    return integrator
end

function get_ode_solution(func, u0, tspan, dt=nothing, param=nothing, alg=Tsit5();
    reltol=1e-12, abstol=1e-12, save_everystep=true)
    """
    Args:
    - f   : the update function
    - u0  : Inital vector of value
    - param : additional param
    - dt  : Time step for integration (when adaptive, it is initial dt)
    - alg : update algorithm (e.g. Euler(), RK4())
    """

    if param == nothing
        prob = ODEProblem(func, u0, tspan)
    else
        prob = ODEProblem(func, u0, tspan, param)
    end

    if dt == nothing
        sol = solve(prob, alg, reltol=reltol, abstol=abstol, adaptive=true, save_everystep=save_everystep)
    else
        sol = solve(prob, alg, dt=dt, reltol=reltol, abstol=abstol, adaptive=false, save_everystep=save_everystep)
    end
    return sol
end
