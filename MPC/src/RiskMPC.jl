using JuMP
using Ipopt
using NamedArrays

using Plots
using JLD2
using Dates

include("structures/nodes.jl")
include("structures/uncertainty.jl")
include("structures/feedback.jl")
include("structures/mpc.jl")


export GasTurbine, Battery, WindTurbine, Demand, WindTurbineSDE, SimpleMarket, Node
export AbstractMPC, StochasticMPC, DeterministicMPC, RiskMeasureMPC
export MPC
export SDE
export Expectation, Entropic, AbstractRiskMeasure, CVaR
export ScenarioTree

include("variables.jl")
include("constraints.jl")
include("cost.jl")
include("ocp.jl")
include("utils/utils.jl")

export set_initial_value, set_actions!, initialize_actions
export build_mpc, set_x₀!, set_profile!, solve!
export build_feedback_struct
export build_scenariotree

include("utils/plotting.jl")
export plot_openloop, plot_closed_loop
export plot_history
