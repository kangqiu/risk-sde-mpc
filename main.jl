using Revise
using Dates
using Plots
using Random
using NamedArrays
using JuMP
using ProgressBars
using Interpolations
using Distributions
using DataInterpolations


include("utils.jl")
include("MPC/src/RiskMPC.jl")
include("Simulator/src/Simulator.jl")

Random.seed!(1234);
Δt = 5; # time step in minutes
T_mpc = 6; # MPC horizon in hours
start = DateTime(2022, 07, 3, 0, 0, 0);
stop = DateTime(2022, 07, 4, 0, 0, 0);
N_sim = Int(Minute(stop-start).value/Δt);  # number of simulation steps
N_mpc = Int(T_mpc*60/Δt); # number of MPC steps
tspan = (0.0, Hour(stop-start).value)
Ω = 9

#### Simulation data
demand = get_simulation_demand(Δt, Day((stop-start)).value+1, 2.7); # demand profile in MW, 3 days
plot(demand, ylims=(0, 10))
spot = get_spot_data(start, stop, T_mpc)
plot(spot["buy"])
price_buy = spot["buy"]
price_sell = spot["sell"]

### here we prep the simulation of wind using SDEs
df_nwp = get_historical_forecast_data(start, stop, T_mpc);
#this calls one solution of the SDE using JuliaSims SDE framework and solver (which we wish to Monte Carlo sample efficiently)
v_wind, p_wind, r_wind, q_wind, p_wind_l, prob = simulate_sde(df_nwp, start, stop, tspan, Δt, N_mpc);

# #you can also simulate multiple trajectories with ensemble problems and analyze the moments of the solution (1000 trajectories)
# @time sol, summ = simulate_ensemble_sde(df_nwp, start, stop, tspan, Δt, 1000);
# plot(summ)


data = Dict(
    "start" => start,
    "stop" => stop,
    "N_sim" => N_sim,
    "T_mpc" => T_mpc,
    "price_buy" => price_buy, 
    "price_sell" => price_sell, 
    "v_wind" => v_wind,
    "power_pot" => p_wind,
    "power_pot_l" => p_wind_l,
    "r_wind" => r_wind,
    "q_wind" => q_wind,
    "demand" => demand,
    "Ω" => Ω);

#weights of the cost function
weights = Dict(
    "balance" => 10.,
    "bat" => 0.1,
    "bat_chargedischarge" => 10.,
    "bat_min" => 0.1,
    "bat_max" => 0.1,
    "bat_terminal" => 0.01,
    "spot" => 0.1,
    "market" => 0.01,
);
nodes = NamedArray(
    [
        # GasTurbine("GT1", 34),
        Battery("B1", 5., 2*60*60., 0.95),
        WindTurbineSDE("WT1", 15),
        Demand("D1"),
        SimpleMarket("M1", 15)
    ],
    [
        # "GT1",
        "B1", 
        "WT1", 
        "D1",
        "M1"
        ]
)

# risk_objective = Dict()

risk = Expectation()
#risk = CVaR(90)
# mpc = DeterministicMPC(Δt, 1:N_mpc, weights)

mpc = RiskMeasureMPC(
    Δt, 
    1:N_mpc,
    weights,
    risk,
    ScenarioTree(Ω)
)



ocp = build_mpc(mpc, nodes)
#mpc = DeterministicMPC(Δt, 1:N_mpc, weights)


feedback = set_initial_value(nodes.array, data);
actions = initialize_actions(collect(nodes));

history = initialize_history()
#unset_silent(ocp)


for k in ProgressBar(1:N_sim)
    # unset_silent(ocp) 
    # k = 10
    set_x₀!(ocp, feedback)

    samples, wind = get_samples(mpc, df_nwp, data, k);

    set_profile!(ocp, nodes["D1"], demand[k:k+N_mpc-1], mpc);
    set_profile!(ocp, nodes["WT1"], samples, mpc);
    set_profile!(ocp, nodes["M1"], price_buy[k:k+N_mpc-1], price_sell[k:k+N_mpc-1], mpc)

    solve!(ocp, true)

    set_actions!(ocp, actions, mpc)
    simulate_system!(data, nodes, actions, feedback, k, mpc)
    update_history!(history, actions, feedback, demand[k:k+N_mpc-1], samples, ocp)

    # plot_openloop(ocp, collect(nodes), mpc)
end

states, inputs = plot_closed_loop(nodes, history, data, "expectation");
