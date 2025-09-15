"""
    set_profile!(ocp::Model, node::Demand, samples::Vector, mpc::DeterministicMPC)

Sets the demand profile for a given node in the optimal control problem (OCP) model using provided sample trajectories.

# Arguments
- `ocp::Model`: The JuMP model representing the optimal control problem.
- `node::Demand`: The demand node for which the profile is set.
- `samples::Vector`: A vector of estimated power for each timestep and percentile (50% is mean since distribution is assumed Gaussian)
- `mpc::DeterministicMPC`: The deterministic MPC object containing the time index set `𝒩`.

# Description
For each time step in the MPC horizon, this function sets the demand parameter in the OCP model for the given node using the median (50th percentile) value from the provided samples.
"""

function set_profile!(ocp::Model, node::WindTurbineSDE, power_samples::Vector, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    for k ∈ 𝒩
        set_parameter_value(ocp[:windpower][node.id, k], power_samples[k])
    end

end

function set_profile!(ocp::Model, node::WindTurbineSDE, power_samples::Vector, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω
    quantiles = mpc.model.quantiles
    
    for k ∈ mpc.𝒩
      for ω in Ω
        set_parameter_value(ocp[:windpower][node.id, ω, k], power_samples[k][quantiles[ω]])
      end
        # set_parameter_value(ocp[:windpower][node.id, 1, k], 0)
        # set_parameter_value(ocp[:windpower][node.id, 2, k], 0.5)
        # set_parameter_value(ocp[:windpower][node.id, 3, k], 1)
    end

end

function set_profile!(ocp::Model, node::Demand, profile::Vector, mpc::AbstractMPC)
    𝒩 = mpc.𝒩

    for k ∈ 𝒩
        set_parameter_value(ocp[:demand][node.id, k], profile[k])
    end
end

function set_profile!(ocp::Model, node::SimpleMarket, price_buy::Vector, price_sell::Vector, mpc::AbstractMPC)
    𝒩 = mpc.𝒩

    for k ∈ 𝒩
        set_parameter_value(ocp[:price_buy][node.id, k], price_buy[k])
        set_parameter_value(ocp[:price_sell][node.id, k], price_sell[k])
    end
end



function set_x₀!(ocp::Model, feedback::NamedArray)
    for ((name,), value) in enamerate(feedback)
        if typeof(value) == BatteryFeedback
            set_parameter_value(ocp[:x₀][name, "soc"], value.soc)
        end
        if typeof(value) == SimpleMarketFeedback
            set_parameter_value(ocp[:x₀][name,"power_buy"], value.power_buy)
            set_parameter_value(ocp[:x₀][name,"power_sell"], value.power_sell)
        end
    end
end




function set_initial_value(nodes::Vector{Node}, data::Dict)
    names = []
    array = []
    for node in nodes
        if typeof(node) == GasTurbine
            push!(names, node.id)
            push!(array, GTFeedback(node.id, 0.0, 0.0, 0.0))
        elseif typeof(node) == Battery
            push!(names, node.id)
            push!(array, BatteryFeedback(node.id, 0.5, 0.0, 0.0))
        elseif typeof(node) == WindTurbineSDE
            push!(names, node.id)
            push!(array, WindFeedback(node.id,
             data["power_pot"][1], 
             data["power_pot"][1], 
             data["v_wind"][1], 
             data["r_wind"][1], 
             data["q_wind"][1]))
        elseif typeof(node) == Demand
            push!(names, node.id)
            push!(array, DemandFeedback(node.id, data["demand"][1]))
        elseif typeof(node) == SimpleMarket
            push!(names, node.id)
            push!(array, SimpleMarketFeedback(node.id, 0.0, 0.0))
        end
    end
    feedback = NamedArray(array, names)
    return feedback 
end


function initialize_actions(nodes::Vector)
    names = []
    array = []

    for node in nodes
        if typeof(node) == GasTurbine
            push!(names, node.id)
            push!(array, GTAction(node.id, 0.0))
        elseif typeof(node) == Battery
            push!(names, node.id)
            push!(array, BatteryAction(node.id, 0.0, 0.0))
        elseif typeof(node) == WindTurbineSDE
            push!(names, node.id)
            push!(array, WindAction(node.id, 0.0))
        elseif typeof(node) == Demand
            push!(names, node.id)
            push!(array, DemandAction(node.id, 0.0))
        elseif typeof(node) == SimpleMarket
            push!(names, node.id)
            push!(array, SimpleMarketAction(node.id, 0.0, 0.0))
        end
    end 
    action = NamedArray(array, names)
    return action
end

function set_actions!(ocp::Model, actions::NamedArray, mpc::AbstractMPC, tolerance::Real=1e-5)

    gt_actions = filter(x -> isa(x, GTAction), actions)
    if length(gt_actions) > 0
        for action in gt_actions
            load = value(ocp[:u][action.id, "load", 1])
            if load > tolerance
                action.load = load
            else
                action.load = 0
            end
        end
    end

    bat_actions = filter(x -> isa(x, BatteryAction), actions)
    if length(bat_actions) > 0
        for action in bat_actions
            charge = value(ocp[:u][action.id, "charge", 1])
            discharge = value(ocp[:u][action.id, "discharge", 1]) 
            if charge > tolerance
                action.charge = charge
            else
                action.charge = 0
            end

            if discharge > tolerance
                action.discharge = discharge
            else
                action.discharge = 0
            end
        end
    end

    wind_actions = filter(x -> isa(x, WindAction), actions)
    if length(wind_actions) > 0
        for action in wind_actions
            curtailment = value(ocp[:u][action.id, "curt", 1])
            if curtailment > tolerance 
                action.curtailment = curtailment
            else
                action.curtailment = 0
            end
        end
    end

    demand_actions = filter(x -> isa(x, DemandAction), actions)
    if length(demand_actions) > 0
        for action in demand_actions
            action.P = value(ocp[:x][action.id, "P", 1])
        end
    end

    simplemarket_actions = filter(x -> isa(x, SimpleMarketAction), actions)
    if length(simplemarket_actions) > 0
        for action in simplemarket_actions
            action.buy = clamp(value(ocp[:u][action.id, "buy", 2]), 0, 1)
            action.sell = clamp(value(ocp[:u][action.id, "sell", 2]), 0, 1)
        end
    end

end

function set_actions!(ocp::Model, actions::NamedArray, mpc::RiskMeasureMPC, tolerance::Real=1e-5)

    gt_actions = filter(x -> isa(x, GTAction), actions)
    if length(gt_actions) > 0
        for action in gt_actions
            load = value(ocp[:u][action.id, "load", 1, 1])
            if load > tolerance
                action.load = load
            else
                action.load = 0
            end
        end
    end

    bat_actions = filter(x -> isa(x, BatteryAction), actions)
    if length(bat_actions) > 0
        for action in bat_actions
            charge = value(ocp[:u][action.id, "charge", 1, 1])
            discharge = value(ocp[:u][action.id, "discharge", 1, 1]) 
            if charge > tolerance
                action.charge = charge
            else
                action.charge = 0
            end

            if discharge > tolerance
                action.discharge = discharge
            else
                action.discharge = 0
            end
        end
    end

    wind_actions = filter(x -> isa(x, WindAction), actions)
    if length(wind_actions) > 0
        for action in wind_actions
            curtailment = value(ocp[:u][action.id, "curt", 1, 1])
            if curtailment > tolerance 
                action.curtailment = curtailment
            else
                action.curtailment = 0
            end
        end
    end

    demand_actions = filter(x -> isa(x, DemandAction), actions)
    if length(demand_actions) > 0
        for action in demand_actions
            action.P = value(ocp[:x][action.id, "P", 1, 1])
        end
    end
    
    simplemarket_actions = filter(x -> isa(x, SimpleMarketAction), actions)
    if length(simplemarket_actions) > 0
        for action in simplemarket_actions
            action.buy = clamp(value(ocp[:u][action.id, "buy", 1, 2]), 0, 1)
            action.sell = clamp(value(ocp[:u][action.id, "sell", 1, 2]), 0, 1)
        end
    end

end

