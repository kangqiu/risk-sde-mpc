using Revise
using Dates
using Plots
using Random
using NamedArrays
using JuMP


using SDESimulator

theta_v_d, theta_R_d, theta_P_d, theta_Q_d, rho_d, mu_d, _, _, xi1_d, xi2_d, xi3_d, sigmaZ, sigma_R_d, sigma_P_d, sigma_Q_d, wp, dwp = params.θ
t = 1
r = 1
vp = 1
q = 1

A1 = [
        ((2*(rho_d*dwp(t)+r)*exp(-v^2/4)*v^2-2*theta_v_d*v^2)/(2*v^2)-4*(rho_d*dwp(t)+ r)*(1-exp(-v^2/4)+theta_v_d*(wp(t)*mu_d-v^2)-sigmaZ^2)/(2*v^2))
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

       Q = [
        sigmaZ 0 0 0; 
        0 sigma_R_d 0 0;
        0 0 sigma_P_d 0;
        0 0 0 sigma_Q_d
    ]

    P = [
        0.0 0.0 0.0 0.0; 
        0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 0.0;
        0.0 0.0 0.0 1.5
    ]
dP = A*P + P*A' + Q*Q'




(4 * (rho_d * dwp(t) + r) * (1 - exp(-v^2 / 4)) + theta_v_d * (4 * wp(t) * mu_d - v^2) - sigmaZ^2) / (2 * v)

function lamperti_ode(dx, x, θ, t)
    theta_v_d, theta_R_d, theta_P_d, theta_Q_d, rho_d, mu_d, _, _, xi1_d, xi2_d, xi3_d, sigmaZ, sigma_R_d, sigma_P_d, sigma_Q_d, wp, dwp = θ
    v, r, vp, q, P = x
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

   A1 = [
        ((2*(rho_d*dwp(t)+r)*exp(-v^2/4)*v^2-2*theta_v_d*v^2)/(2*v^2)-4*(rho_d*dwp(t)+ r)*(1-exp(-v^2/4)+theta_v_d*(wp(t)*mu_d-v^2)-sigmaZ^2)/(2*v^2))
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

       Q = [
        sigmaZ 0 0 0; 
        0 sigma_R_d 0 0;
        0 0 sigma_P_d 0;
        0 0 0 sigma_Q_d
    ]

    dP = A*P + P*A' + Q*Q'
    dx[5] = dP[1,1]
    dx[6] = dP[1,2]
    dx[7] = dP[1,3]
    dx[8] = dP[1,4]
    dx[9] = dP[2,1]
    dx[10] = dP[2,2]
    dx[11] = dP[2,3]
    dx[12] = dP[2,4]
    dx[13] = dP[3,1]
    dx[14] = dP[3,2]
    dx[15] = dP[3,3]
    dx[16] = dP[3,4]
    dx[17] = dP[4,1]
    dx[18] = dP[4,2]
    dx[19] = dP[4,3]
    dx[20] = dP[4,4]
end

using OrdinaryDiffEq, Plots
gr()



#Half-life of Carbon-14 is 5,730 years.
t½ = 5.730

#Setup
u₀ = 1.0
tspan = (0.0, 30.0) #in hours

#Define the problem
radioactivedecay(u, p, t) = -log(2) / t½ * u

#Pass to solver
prob = ODEProblem(radioactivedecay, u₀, tspan)
sol = solve(prob, Tsit5())

#Plot
plot(sol, linewidth = 2, title = "Carbon-14 half-life",
    xaxis = "Time in thousands of years", yaxis = "Ratio left",
    label = "Numerical Solution")
plot!(sol.t, t -> 2^(-t / t½), lw = 3, ls = :dash, label = "Analytical Solution")

Δt = 5 # time step in minutes
T_mpc = 24 *(60 ) # MPC horizon in minutes
start = DateTime(2022, 07, 3, 12, 0, 0)
stop = DateTime(2022, 07, 5, 12, 0, 0)
N_sim = Int(Minute(stop-start).value/Δt)  # number of simulation steps
N_mpc = Int(T_mpc/Δt)  # number of MPC steps
tspan = (0.0, Hour(stop-start).value)
tspan = (0.0, 24.0) # in hours
forecast = get_mpc_nwp(Δt, start, stop, tspan, Int(T_mpc));
sim_param = SimulationParameters(Δt, start, stop)
mpc_time = MPCTimestruct(Δt, Int(T_mpc/Δt))

_, df_nwp = get_nwp(
    "/Users/kangqiu/Documents/PhD/Github/SDEMonteCarlo/data/winddata/data_forecast.jld2", 
    sim_param, mpc_time)

sim_nwp, sim_dnwp, sim_x0 = get_simulation_nwp(df_nwp, sim_param, tspan)

params = SDEParameters(sim_nwp, sim_dnwp)

theta_v_d, theta_R_d, theta_P_d, theta_Q_d, rho_d, mu_d, _, _, xi1_d, xi2_d, xi3_d, sigmaZ, sigma_R_d, sigma_P_d, sigma_Q_d, wp, dwp = params.θ

prob = ODEProblem(lamperti_ode, x₀, tspan, params.θ)

sol = solve(prob, Tsit5())
plot(sim_nwp.(sol.t))


v_l = [s[1] for s in sol.u]
plot(v_l, label = "Wind speed in Lamperti domain", xlabel = "Time (hours)", ylabel = "Wind speed (m/s)", title = "Wind Speed in Lamperti Domain")
pv = [s[3] for s in sol.u]
plot(pv, label = "Wind power in Lamperti domain", xlabel = "Time (hours)", ylabel = "Wind power", title = "Wind Power in Lamperti Domain")

P_pred = [s[5:20] for s in sol.u]
P_pred = reshape.(P_pred, 4, 4)
σ_p = [P_pred[i][3,3] for i in 1:length(P_pred)]
plot(σ_p, label = "Predicted variance of wind power", xlabel = "Time (hours)", ylabel = "Variance", title = "Predicted Variance of Wind Power")