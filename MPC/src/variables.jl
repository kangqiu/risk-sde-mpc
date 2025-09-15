function define_variables!(ocp::Model, nodes::NamedArray, mpc::DeterministicMPC)
    #collect indices
    𝒩 = mpc.𝒩

    name = names(nodes)[1]
    x_values = []
    u_values = []

    for node in nodes.array
        x_ind, u_ind = define_variable_indeces!(node)
        push!(x_values, x_ind)
        push!(u_values, u_ind)
    end
    x_dict = Dict(zip(name, x_values))
    u_dict = Dict(zip(name, u_values))

    @variable(ocp, x[i=name, x_dict[i], 𝒩])
    @variable(ocp, u[i=name, u_dict[i], 𝒩])
    @variable(ocp, x₀[i=name, x_dict[i]] in Parameter(0.0))  # initial values

    #define slack variables
    σ_values = []
    σ_name =[]
    for node in nodes.array
        if typeof(node) <: Battery
            #so far only Battery nodes have slack variables
            push!(σ_values, define_slack_variable_indeces!(node))
            push!(σ_name, node.id)
        end
    end

    σ_dict = Dict(zip(σ_name, σ_values))

    @variable(ocp, σ[i=σ_name, σ_dict[i], 𝒩])
    @variable(ocp, balance[𝒩])  # balance variable for energy balance constraints
end

function define_variables!(ocp::Model, nodes::NamedArray, mpc::RiskMeasureMPC)
    #collect indices
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω

    name = names(nodes)[1]
    x_values = []
    u_values = []

    for node in nodes.array
        x_ind, u_ind = define_variable_indeces!(node)
        push!(x_values, x_ind)
        push!(u_values, u_ind)
    end
    x_dict = Dict(zip(name, x_values))
    u_dict = Dict(zip(name, u_values))

    @variable(ocp, x[i=name, x_dict[i], Ω, 𝒩])
    @variable(ocp, u[i=name, u_dict[i], Ω, 𝒩])
    @variable(ocp, x₀[i=name, x_dict[i]] in Parameter(0.0))  # initial values

    #define slack variables
    σ_values = []
    σ_name =[]
    for node in nodes.array
        if typeof(node) <: Battery
            #so far only Battery nodes have slack variables
            push!(σ_values, define_slack_variable_indeces!(node))
            push!(σ_name, node.id)
        end
    end

    σ_dict = Dict(zip(σ_name, σ_values))

    @variable(ocp, σ[i=σ_name, σ_dict[i], Ω, 𝒩])
    @variable(ocp, balance[Ω, 𝒩])  # balance variable for energy balance constraints
end


function define_variable_indeces!(node::GasTurbine)
   return (x = ["η", "CO₂", "power"], u = ["load"])
end

function define_variable_indeces!(node::Battery)
   return (x = ["soc", "Pₒᵤₜ", "Pᵢₙ"], u = ["charge", "discharge"])
end

function define_variable_indeces!(node::WindTurbineSDE)
   return (x = ["power"], u = ["curt"])
end

function define_variable_indeces!(node::Demand)
   return (x = ["P"], u = ["ΔP"])
end

function define_variable_indeces!(node::SimpleMarket)
   return (x = ["power_buy", "power_sell"], u = ["buy", "sell"])
end

function define_parameters!(ocp::Model, nodes::NamedArray, mpc::AbstractMPC)

    WTSDE_nodes = filter(x -> typeof(x)==WindTurbineSDE, nodes.array)
    WTSDE_nodes = map(x -> x::WindTurbineSDE, WTSDE_nodes) #map to a narrower array type for multiple dispatch
    if length(WTSDE_nodes) > 0 
        define_parameters!(ocp, WTSDE_nodes, mpc) 
    end

    demand_nodes = filter(x -> typeof(x)==Demand, nodes.array)
    demand_nodes = map(x -> x::Demand, demand_nodes) #map to a narrower array type for multiple dispatch
    if length(demand_nodes) > 0
        define_parameters!(ocp, demand_nodes, mpc)  
    end

    market_nodes = filter(x -> typeof(x)==SimpleMarket, nodes.array)
    market_nodes = map(x -> x::SimpleMarket, market_nodes) #map to a narrower array type for multiple dispatch
    if length(market_nodes) > 0
        define_parameters!(ocp, market_nodes, mpc)  
    end

end

function define_parameters!(ocp::Model, nodes::Array{WindTurbineSDE}, mpc::AbstractMPC)
    𝒩 = mpc.𝒩
    ids = []
    for n in nodes
        push!(ids, n.id)
    end
    @variable(ocp, windpower[ids, 𝒩] in Parameter(0.0))

end

function define_parameters!(ocp::Model, nodes::Array{WindTurbineSDE}, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω
    ids = []
    for n in nodes
        push!(ids, n.id)
    end
    @variable(ocp, windpower[ids, Ω, 𝒩] in Parameter(0.0))

end

function define_parameters!(ocp::Model, nodes::Array{Demand}, mpc::AbstractMPC)
    𝒩 = mpc.𝒩
    ids = []
    for n in nodes
        push!(ids, n.id)
    end
    @variable(ocp, demand[ids, 𝒩] in Parameter(0.0))
end

function define_parameters!(ocp::Model, nodes::Array{SimpleMarket}, mpc::AbstractMPC)
    𝒩 = mpc.𝒩
    ids = []
    for n in nodes
        push!(ids, n.id)
    end
    @variable(ocp, price_buy[ids, 𝒩] in Parameter(0.0))
    @variable(ocp, price_sell[ids, 𝒩] in Parameter(0.0))
end


function define_slack_variable_indeces!(node::Battery)
   return (σ = ["socₘᵢₙ", "socₘₐₓ"])
end
