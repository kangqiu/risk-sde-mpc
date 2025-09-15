function drift!(dx, x, θ, t)
    v, r, P_wtg, q = x
    theta_v_d, theta_R_d, theta_P_d, theta_Q_d, rho_d, mu_d, _, _, xi1_d, xi2_d, xi3_d, sigmaZ, sigma_R_d, sigma_P_d, sigma_Q_d, wp, dwp = θ

    dx[1] = dv = (4 * (rho_d * dwp(t) + r) * (1 - exp(-v^2 / 4)) + theta_v_d * (4 * wp(t) * mu_d - v^2) - sigmaZ^2) / (2 * v)
    dx[2] = dr = -theta_R_d * r
    dx[3] = dvp = (2 * theta_P_d * (xi3_d / (1 + exp(-xi1_d * (v^ 2 / 4 - xi2_d + q))) - 0.5 * (1 + sin(P_wtg)))) / cos(P_wtg)
    dx[4] = dq = -theta_Q_d * q
end

function diffusion!(dx, x, θ, t)
    theta_v_d, theta_R_d, theta_P_d, theta_Q_d, rho_d, mu_d, _, _, xi1_d, xi2_d, xi3_d, sigmaZ, sigma_R_d, sigma_P_d, sigma_Q_d, wp, dwp = θ
    dx[1] = sigmaZ
    dx[2] = sigma_R_d
    dx[3] = sigma_P_d
    dx[4] = sigma_Q_d
end

function lamperti_covariance_matrix(x, θ, t)
     theta_v_d, theta_R_d, theta_P_d, theta_Q_d, rho_d, mu_d, _, _, xi1_d, xi2_d, xi3_d, sigmaZ, sigma_R_d, sigma_P_d, sigma_Q_d, wp, dwp = θ
     v = x[1]
    r = x[2]
    vp = x[3]
    q = x[4]
    
    A1 = [
        ((2*(rho_d*dwp(t)+r)*exp(-v^2/4)*v^2 - 2 * theta_v_d * v^2)/(2*v^2)-4*(rho_d*dwp(t)+ r)*(1-exp(-v^2/4)+theta_v_d*(wp(t)*mu_d-v^2)-sigmaZ^2)/(2*v^2))
        ((4 - 4 * exp(-v^2 / 4)) / (2 * v))
        0
        0]

    A2 = [0 
        -theta_R_d 
        0 
        0]

    A3 = [theta_P_d * xi3_d * xi1_d * v * exp(-xi1_d * (v^2 / 4 + q - xi2_d)) / (
                                (1 + exp(-xi1_d * (v^2 / 4 + q - xi2_d)))^2 * cos(vp))
            0
            -theta_P_d + 2 * theta_P_d * (
                                    xi3_d / (1 + exp(-xi1_d * (v^2 / 4 + q - xi2_d))) - 0.5 - 0.5 * sin(vp)) * sin(
                            vp) / cos(vp)^2
            2 * theta_P_d * xi3_d * xi1_d * exp(-xi1_d * (v^2 / 4 + q - xi2_d)) / (
                                    (1 + exp(-xi1_d * (v^2 / 4 + q - xi2_d)))^2 * cos(vp))]
        
    A4 =  [0, 0, 0, -theta_Q_d]

    return hcat(A1, A2, A3, A4)
end



"""
    simulate_sde(x₀, tspan, wp, dwp)

Simulates the SDE system for wind power using the given initial state, time span, and wind profile functions.

Arguments:
    x₀::Vector{<:Real}: Initial state vector in natural domain
    tspan::Tuple{<:Real,<:Real}: Time span for simulation
    wp::Any: Wind profile function
    dwp::Any: Derivative of wind profile function

Returns:
    prob: SDEProblem instance
    sol: Solution object
"""

function simulate_sde(x₀::Vector{<:Real}, tspan::Tuple{<:Real,<:Real}, wp::Any,  dwp::Any)
    params = SDEParameters(wp, dwp)
    prob = SDEProblem(drift!, diffusion!, x₀, tspan, params.θ)
    sol = StochasticDiffEq.solve(prob, SOSRI())
    return prob, sol
end


"""
    natural_to_lamperti_transform(x_n)

Transforms the state from the natural domain to the Lamperti domain.

Arguments:
    x_n::Vector: State in natural domain

Returns:
    x_l::Vector: State in Lamperti domain
"""
function natural_to_lamperti_transform(x_n)
    v_wind_l = 2*sqrt(x_n[1])

    x_l = [
        v_wind_l,
        x_n[2],
        # x_n[3],
        asin(x_n[3]/2 -1),
        x_n[4]
        ]

    return x_l
end

function natural_to_lamperti_transform_power(x_n)

    x_l = asin(x_n/2 -1)

    return x_l
end

function lamperti_to_natural_transform_wind(x_l::Real)
    return(x_l^2/4)
end


function  lamperti_to_natural_transform_power(x_l::Real, data_l::Real, data_n::Real)
    Δsample_l = x_l + data_l

    Δpower_n =  (0.5*(1+sin((x_l+data_l)/2))-(0.5*(1+sin(data_l))))*2

    power_n = data_n + Δpower_n
    power_n = clamp(power_n, 0, 1)
    # if α >= 0.5
    #     power_n = data_n + Δpower_n
    # else
    #     power_n = data_n - Δpower_n
    # end

    return(power_n)
end
"""
    lamperti_to_natural_transform(x_l)

Transforms the state from the Lamperti domain back to the natural domain.

Arguments:
    x_l::Vector: State in Lamperti domain

Returns:
    x_n::Vector: State in natural domain
"""
function lamperti_to_natural_transform(x_l::Vector)
    x_n =  [
    x_l[1]^2/4,
    x_l[2],
    0.5*(1+sin(x_l[3])), 
    x_l[4]
    ]
    return x_n
end

function lamperti_to_natural_transform_power(power_l::Real)
    # return 2*(1+sin(power_l))
    return 0.5*(1+sin(power_l))
end

"""
    simulate_sde(
        df_nwp::DataFrame,
        start::DateTime,
        stop::DateTime,
        tspan::Tuple,
        Δt::Int,
        plots=true
    )

Simulates the stochastic differential equation (SDE) system for wind power using historical NWP (Numerical Weather Prediction) data 
and returns interpolated wind speed and power output on a regular time grid.

# Arguments
- `df_nwp::DataFrame`: DataFrame containing NWP data with "forecast" and "dnwp" columns.
- `start::DateTime`: Simulation start time.
- `stop::DateTime`: Simulation stop time.
- `tspan::Tuple`: Tuple specifying the time span for the SDE solver.
- `Δt::Int`: Time step in minutes for the output time grid.
- `plots::Bool=true`: If true, plots wind speed and power output over time.

# Returns
- `v::Vector`: Interpolated wind speed [m/s] over the simulation period in natural domain.
- `p::Vector`: Interpolated normalized power output [⋅] over the simulation period in natural domain.

# Description
This function creates interpolation objects for wind forecast and its increments, initializes the SDE state, runs the SDE simulation, transforms the solution to the natural domain, and optionally plots wind speed and power output. The results are interpolated onto a regular time grid defined by `Δt`.
"""
function simulate_sde(df_nwp::DataFrame,start::DateTime, stop::DateTime, tspan::Tuple, Δt::Int, N_mpc::Int, plots=true)
    
    itp_nwp, itp_dnwp = get_interp_object_best_nwp(df_nwp, start, stop,  Δt,  N_mpc)

    v_wind = itp_nwp(0) # initial wind speed in m/s
    x₀ = [
        v_wind,
        0.0,
        wind_curve(v_wind),
        0.0
    ]
    
    x₀_l = natural_to_lamperti_transform(x₀)

    # here comes the weird stuff the sde takes the natural domain power as input, 
    # since also the transform is not an exact inverse

    x₀ = [
        x₀_l[1],
        x₀_l[2],
        x₀[3],
        x₀_l[4],
    ]

    tspan_adjusted = (tspan[1]-0.9, tspan[2]+Hour(Minute(Δt)*N_mpc).value)
    prob, sol = simulate_sde(x₀, tspan_adjusted, itp_nwp, itp_dnwp)

    # Transform solution to natural domain 
    x_n = [lamperti_to_natural_transform(sol.u[i]) for i in 1:length(sol.u)]
    v_n = [x[1] for x in x_n]
    p_n = [x[3] for x in x_n]
    r_n = [x[2] for x in x_n]
    q_n = [x[4] for x in x_n]

    x_l = [sol.u[i] for i in 1:length(sol.u)]
    p_l = [x[3] for x in x_l]

    #again the first power value is already in the natural domain
    p_n[1] = sol.u[1][3]

    itp_v = DataInterpolations.LinearInterpolation(v_n, sol.t)
    itp_p = DataInterpolations.LinearInterpolation(p_n, sol.t)
    itp_r = DataInterpolations.LinearInterpolation(r_n, sol.t)
    itp_q = DataInterpolations.LinearInterpolation(q_n, sol.t)
    itp_p_l = DataInterpolations.LinearInterpolation(p_l, sol.t)



    timegrid = tspan[1]:Δt/60:tspan[2]+Hour(Minute(Δt)*N_mpc).value
    v = itp_v.(timegrid)
    p = itp_p.(timegrid)
    r = itp_r.(timegrid)
    q = itp_q.(timegrid)
    p_l = itp_p_l.(timegrid)
    
    if plots
        p1 = plot(timegrid, v, label = "wind speed [m/s]", xlabel = "time [h]")
        p2 = plot(timegrid, p, label = "power output [⋅]",  xlabel = "time [h]")

        display(plot(p1, p2, layout =(2,1), title="wind speed and power"))
    end


    return v, p, r, q, p_l, prob
end

function simulate_ensemble_sde(df_nwp::DataFrame, x₀::Array, start::DateTime, stop::DateTime, tspan::Tuple, Δt::Int, trajectories::Int)
    itp_nwp, itp_dnwp = get_interp_object_best_nwp(df_nwp, start, stop,  Δt,  0)
    
    x₀_l = natural_to_lamperti_transform(x₀)

    # here comes the weird stuff the sde takes the natural domain power as input, 
    # since also the transform is not an exact inverse

    x₀ = [
        x₀_l[1],
        x₀_l[2],
        x₀[3],
        x₀_l[4],
    ]

    tspan = (tspan[1], tspan[2])
    params = SDEParameters(itp_nwp, itp_dnwp)
    prob = SDEProblem(drift!, diffusion!, x₀, tspan, params.θ)

    ensembleprob = EnsembleProblem(prob)
    sol = StochasticDiffEq.solve(ensembleprob, SOSRI(), trajectories=trajectories)
    summ = EnsembleSummary(sol, tspan[1]:Δt/60:tspan[2])
    # display(plot(summ))
    return sol
end

function lamperti_ode(dx, x, θ, t)
    (theta_v_d, 
    theta_R_d,
    theta_P_d,
    theta_Q_d, 
    rho_d,
    mu_d, _, _, xi1_d, xi2_d, xi3_d, sigmaZ, sigma_R_d, sigma_P_d, sigma_Q_d, wp, dwp) = θ
    v = x[1]
    r = x[2]
    vp = x[3]
    q = x[4]
    P = reshape(x[5:end], 4, 4)
    # Lamperti ODE system
    dx[1] = dv = (4 * (rho_d * dwp(t) + r) * (1 - exp(-v^2 / 4)) + theta_v_d * (4 * wp(t) * mu_d - v^2) - sigmaZ^2) / (2 * v)
    dx[2] = dr = -theta_R_d * r
    dx[3] = dvp = (2 * theta_P_d * (xi3_d / (1 + exp(-xi1_d * (v^ 2 / 4 - xi2_d + q))) - 0.5 * (1 + sin(vp)))) / cos(vp)
    dx[4] = dq = -theta_Q_d * q

    Q = [
        sigmaZ 0 0 0; 
        0 sigma_R_d 0 0;
        0 0 sigma_P_d 0;
        0 0 0 sigma_Q_d
    ]

    # A = lamperti_covariance_matrix(x, θ, t)
    A1 = [
        ((2*(rho_d*dwp(t)+r)*exp(-v^2/4)*v^2 - 2 * theta_v_d * v^2)/(2*v^2)-4*(rho_d*dwp(t)+ r)*(1-exp(-v^2/4)+theta_v_d*(wp(t)*mu_d-v^2)-sigmaZ^2)/(2*v^2))
        ((4 - 4 * exp(-v^2 / 4)) / (2 * v))
        0
        0]

    A2 = [0 
        -theta_R_d 
        0 
        0]

    A3 = [theta_P_d * xi3_d * xi1_d * v * exp(-xi1_d * (v^2 / 4 + q - xi2_d)) / (
                                (1 + exp(-xi1_d * (v^2 / 4 + q - xi2_d)))^2 * cos(vp))
            0
            -theta_P_d + 2 * theta_P_d * (
                                    xi3_d / (1 + exp(-xi1_d * (v^2 / 4 + q - xi2_d))) - 0.5 - 0.5 * sin(vp)) * sin(
                            vp) / cos(vp)^2
            2 * theta_P_d * xi3_d * xi1_d * exp(-xi1_d * (v^2 / 4 + q - xi2_d)) / (
                                    (1 + exp(-xi1_d * (v^2 / 4 + q - xi2_d)))^2 * cos(vp))]
        
    A4 =  [0, 0, 0, -theta_Q_d]

    A = hcat(A1, A2, A3, A4)

    dx[5:end] = dP = reshape(A*P + P*A'+ Q*Q', 16, 1)

end

function simulate_ode(x₀::Vector{<:Real}, tspan::Tuple{<:Real,<:Real}, wp::Any,  dwp::Any)
    params = SDEParameters(wp, dwp);
    prob = ODEProblem(lamperti_ode, x₀, tspan, params.θ);
    sol = OrdinaryDiffEq.solve(prob, Tsit5());
    return prob, sol
end

function simulate_ode(
    T_mpc::Real,
    data::Dict,
    itp_nwp::DataInterpolations.ConstantInterpolation,
    itp_dnwp::DataInterpolations.ConstantInterpolation,
    Δt::Int,
    k::Int
    )

    u₀ = [
        data["v_wind"][k],
        data["r_wind"][k],
        data["power_pot"][k],
        data["q_wind"][k]
    ];

    u₀_l = natural_to_lamperti_transform(u₀)

    # ODE same as SDE takes natural power instead of lamperti power
     u₀ = [
        u₀_l[1],
        u₀_l[2],
        # u₀[3],
        data["power_pot_l"][k],
        u₀_l[4],
    ];
    
    σ₀ = repeat([0.0], 16);
    x₀ = vcat(u₀, σ₀);

    tspan = (0, T_mpc);

    prob, sol = simulate_ode(x₀, tspan, itp_nwp, itp_dnwp);

    # pick out the mean and covariance of power output in lamperti domain   
    μ_v_l = [x[1] for x in sol.u];
    μ_p_l = [x[3] for x in sol.u];
    σ_v_l = [sqrt(x[5]) for x in sol.u]; 
    σ_p_l = [sqrt(x[15])*10 for x in sol.u]; 


    itp_μ = DataInterpolations.LinearInterpolation(μ_p_l, sol.t);
    itp_μv = DataInterpolations.LinearInterpolation(μ_v_l, sol.t);
    itp_σ = DataInterpolations.LinearInterpolation(σ_p_l, sol.t);
    itp_μσ = DataInterpolations.LinearInterpolation(μ_v_l, sol.t);

    timegrid = 0:Δt/60:T_mpc;
    μ = itp_μ.(timegrid);
    σ = itp_σ.(timegrid);
    μv = itp_μv.(timegrid)
    σv = itp_μσ.(timegrid)
    
    return μ, σ, μv, σv

end 
