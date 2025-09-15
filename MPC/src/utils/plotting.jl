function plot_energy_balance_openloop(ocp::Model, nodes::Vector, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    cumulative_output = zeros(last(𝒩))
    cumulative_input = zeros(last(𝒩))
    p = plot(title = "Open loop energybalance")
    for n in nodes 
        if typeof(n) <: GasTurbine
            gt_power = [cumulative_output[k] + value(ocp[:x][n.id, "power", k]) for k in 𝒩]
            plot!(𝒩, gt_power , fillrange = cumulative_output, palette = :tab10,label = n.id)
            cumulative_output = gt_power 
        elseif typeof(n) <: WindTurbineSDE
            wt_power = [cumulative_output[k] + value(ocp[:x][n.id, "power", k]) for k in 𝒩]
            plot!(𝒩, wt_power, fillrange = cumulative_output, palette = :tab10,label = n.id)
            cumulative_output = wt_power 
        elseif typeof(n) <: Demand
            demand_power = [cumulative_input[k] + value(ocp[:x][n.id, "P", k]) for k in 𝒩]
            plot!(𝒩, -demand_power, fillrange = -cumulative_input, palette = :tab10,label = n.id)
            cumulative_input = demand_power
        elseif typeof(n) <: Battery
            bat_power_in = [cumulative_input[k] + value(ocp[:x][n.id, "Pᵢₙ", k]) for k in 𝒩]
            bat_power_out = [cumulative_output[k] + value(ocp[:x][n.id, "Pₒᵤₜ", k]) for k in 𝒩]

            plot!(𝒩, bat_power_out, fillrange = cumulative_output, palette = :tab10,label = n.id)
            cumulative_output = bat_power_out
            plot!(𝒩, -bat_power_in, fillrange = -cumulative_input, palette = :tab10,label = n.id)
            cumulative_input = bat_power_in 

        elseif typeof(n) <: SimpleMarket
            market_sell = [cumulative_input[k] + value(ocp[:x][n.id, "power_sell", k]) for k in 𝒩]
            market_buy = [cumulative_output[k] + value(ocp[:x][n.id, "power_buy", k]) for k in 𝒩]

            plot!(𝒩, market_sell, fillrange = cumulative_output, palette = :tab10,label = "$n.id sell")
            cumulative_output = market_sell
            plot!(𝒩, -market_buy, fillrange = -cumulative_input, palette = :tab10,label = "$n.id buy")
            cumulative_input = market_buy
        end
    end

    display(p)
end

function plot_energy_balance_openloop(ocp::Model, nodes::Vector, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω

    plots = []
    for ω ∈ Ω
        cumulative_output = zeros(last(𝒩))
        cumulative_input = zeros(last(𝒩))

        p = plot(legend=:outertopright)
        for n in nodes 
            if typeof(n) <: GasTurbine
                gt_power = [cumulative_output[k] + value(ocp[:x][n.id, "power", ω, k]) for k in 𝒩]
                plot!(𝒩, gt_power , fillrange = cumulative_output, fillalpha=0.2, palette = :tab10,label = n.id)
                cumulative_output = gt_power 
            elseif typeof(n) <: WindTurbineSDE
                wt_power = [cumulative_output[k] + value(ocp[:x][n.id, "power", ω, k]) for k in 𝒩]
                plot!(𝒩, wt_power, fillrange = cumulative_output, fillalpha=0.2, palette = :tab10,label = n.id)
                cumulative_output = wt_power 
            elseif typeof(n) <: Demand
                demand_power = [cumulative_input[k] + value(ocp[:x][n.id, "P", ω, k]) for k in 𝒩]
                plot!(𝒩, -demand_power, fillrange = -cumulative_input, fillalpha=0.2, palette = :tab10,label = n.id)
                cumulative_input = demand_power
            elseif typeof(n) <: Battery
                bat_power_in = [cumulative_input[k] + value(ocp[:x][n.id, "Pᵢₙ", ω, k]) for k in 𝒩]
                bat_power_out = [cumulative_output[k] + value(ocp[:x][n.id, "Pₒᵤₜ", ω, k]) for k in 𝒩]

                plot!(𝒩, bat_power_out, fillrange = cumulative_output, fillalpha=0.2, palette = :tab10,label = n.id)
                cumulative_output = bat_power_out
                plot!(𝒩, -bat_power_in, fillrange = -cumulative_input, fillalpha=0.2, palette = :tab10,label = n.id)
                cumulative_input = bat_power_in 
            elseif typeof(n) <: SimpleMarket
                market_sell = [cumulative_input[k] + value(ocp[:x][n.id, "power_sell", ω, k]) for k in 𝒩]
                market_buy = [cumulative_output[k] + value(ocp[:x][n.id, "power_buy", ω, k]) for k in 𝒩]

                plot!(𝒩, market_sell, fillrange = cumulative_output, palette = :tab10,label = n.id)
                cumulative_output = market_sell
                plot!(𝒩, -market_buy, fillrange = -cumulative_input, palette = :tab10,label = n.id)
                cumulative_input = market_buy
            end
        end
        # push!(plots, p)
        display(plot(p, title="openloop energy balance scenario $ω"))

    end
    
    # display(plot(plots..., layout = (length(plots), 1), title="openloop energy balance"))

end

function plot_openloop(ocp::Model, n::Battery, mpc::DeterministicMPC)
    #battery soc
    𝒩 = mpc.𝒩
    soc = [value(ocp[:x][n.id, "soc", k]) for k in 𝒩]
    charge = [value(ocp[:u][n.id, "charge", k]) for k in 𝒩]
    discharge = [value(ocp[:u][n.id, "discharge", k]) for k in 𝒩]

    p1 = plot(𝒩, soc, palette = :tab10, label = "SOC")
    p2 = plot(𝒩, charge, palette = :tab10, label = "charge")
    p2 = plot!(𝒩, discharge, palette = :tab10, label = "discharge")
    display(plot(p1, p2, layout =(2,1), title=n.id))
end

function plot_openloop(ocp::Model, n::Battery, mpc::RiskMeasureMPC)
    #battery soc
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω

    plots = []
    for ω ∈ Ω
        soc = [value(ocp[:x][n.id, "soc", ω, k]) for k in 𝒩]
        charge = [value(ocp[:u][n.id, "charge", ω,  k]) for k in 𝒩]
        discharge = [value(ocp[:u][n.id, "discharge",  ω, k]) for k in 𝒩]

        p1 = plot(𝒩, soc, palette = :tab10, label = "SOC", legend=:outertopright)
        p2 = plot(𝒩, charge, palette = :tab10, label = "charge")
        p2 = plot!(𝒩, discharge, palette = :tab10, label = "discharge", legend=:outertopright)

        display(plot(p1, p2, layout =(2,1), title="scenario $ω" ))
    end
end

function plot_openloop(ocp::Model, n::GasTurbine, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    power = [value(ocp[:x][n.id, "power", k]) for k in 𝒩]
    efficiency = [value(ocp[:x][n.id, "η", k]) for k in 𝒩]
    load = [value(ocp[:u][n.id, "load", k]) for k in 𝒩]

    p1 = plot(𝒩, power, palette = :tab10, label = "power",legend=:outertopright)
    p2 = plot(𝒩, efficiency, palette = :tab10, label = "efficiency",legend=:outertopright)
    p2 = plot!(𝒩, load, palette = :tab10, label = "load")

    display(plot(p1, p2, layout =(2,1), title = n.id))
end

function plot_openloop(ocp::Model, n::GasTurbine, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω

    plots = []
    for ω ∈ Ω
        power = [value(ocp[:x][n.id, "power", ω, k]) for k in 𝒩]
        efficiency = [value(ocp[:x][n.id, "η", ω, k]) for k in 𝒩]
        load = [value(ocp[:u][n.id, "load",  ω, k]) for k in 𝒩]

        p1 = plot(𝒩, power, palette = :tab10, label = "power", legend=:outertopright)
        p2 = plot(𝒩, efficiency, palette = :tab10, label = "efficiency", legend=:outertopright)
        p2 = plot!(𝒩, load, palette = :tab10, label = "load")

        display(plot(p1, p2, layout =(2,1), title = "scenario $ω"))
    end
end


function plot_openloop(ocp::Model, n::WindTurbineSDE, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    power = [value(ocp[:x][n.id, "power", k]) for k in 𝒩]
    power_pot = [n.Pₘₐₓ * value( ocp[:windpower]["WT1", k]) for k in 𝒩]
    curtailed = [n.Pₘₐₓ * value(ocp[:u][n.id, "curt", k]) for k in 𝒩]

    p1 = plot(𝒩, power, palette = :tab10,  fillrange = zeros(last(𝒩)), fillalpha = 0.15, label = "wind power", legend=:outertopright)
    p1 = plot!(𝒩, power_pot, palette = :tab10,  fillrange = power,  fillalpha = 0.15, label = "wind power potential", legend=:outertopright)
    p1 = plot!(𝒩, -curtailed , fillrange = zeros(last(𝒩)), fillalpha = 0.15, palette = :tab10, label = "wind curtailed")
    display(p1)
end

function plot_openloop(ocp::Model, n::WindTurbineSDE, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω

    plots = []

    for ω ∈ Ω
        power = [value(ocp[:x][n.id, "power", ω, k]) for k in 𝒩]
        power_pot = [n.Pₘₐₓ * value( ocp[:windpower]["WT1", ω, k]) for k in 𝒩]
        curtailed = [n.Pₘₐₓ * value(ocp[:u][n.id, "curt",  ω, k]) for k in 𝒩]

        p1 = plot(𝒩, power, palette = :tab10,  fillrange = zeros(last(𝒩)), fillalpha = 0.15, label = "wind power", legend=:outertopright)
        p1 = plot!(𝒩, power_pot, palette = :tab10,  fillrange = power,  fillalpha = 0.15, label = "wind power potential", legend=:outertopright)
        p1 = plot!(𝒩, -curtailed , fillrange = zeros(last(𝒩)), fillalpha = 0.15, palette = :tab10, label = "wind curtailed")
        # push!(plots, p1)
        display(plot(p1, title = "$n.id scenario $ω"))
    end

    # display(plot(plots..., layout =(length(plots),1), title = "$n.id"))
end

function plot_openloop(ocp::Model, n::Demand, mpc::AbstractMPC)
    𝒩 = mpc.𝒩
end


function plot_openloop(ocp::Model, n::SimpleMarket, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    power_buy = [value(ocp[:x][n.id, "power_buy", k]) for k in 𝒩]
    power_sell = [value(ocp[:x][n.id, "power_sell", k]) for k in 𝒩]
    price_buy = [value(ocp[:price_buy][n.id, k]) for k in 𝒩]
    price_sell = [value(ocp[:price_sell][n.id, k]) for k in 𝒩]

    p1 = plot(𝒩, -power_buy, palette = :tab10,  fillrange = zeros(last(𝒩)), fillalpha = 0.15, label = "power bought")
    p1 = plot!(𝒩, power_sell, palette = :tab10,  fillrange = zeros(last(𝒩)), fillalpha = 0.15, label = "power sold")
    p2 = plot(𝒩, -price_buy, palette = :tab10, label = "price buy")
    p2 = plot!(𝒩, price_sell, palette = :tab10, label = "price sell")
    display(plot(p1, p2, layout =(2,1), title=n.id))
end
function plot_openloop(ocp::Model, n::SimpleMarket, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω
    plots = []
    for ω ∈ Ω
        power_buy = [value(ocp[:x][n.id, "power_buy", ω, k]) for k in 𝒩]
        power_sell = [value(ocp[:x][n.id, "power_sell", ω, k]) for k in 𝒩]
        price_buy = [value(ocp[:price_buy][n.id, k]) for k in 𝒩]
        price_sell = [value(ocp[:price_sell][n.id, k]) for k in 𝒩]

        p1 = plot(𝒩, -power_buy, palette = :tab10,  fillrange = zeros(last(𝒩)), fillalpha = 0.15, label = "power bought")
        p1 = plot!(𝒩, power_sell, palette = :tab10,  fillrange = zeros(last(𝒩)), fillalpha = 0.15, label = "power sold")
        p2 = plot(𝒩, -price_buy, palette = :tab10, label = "price buy")
        p2 = plot!(𝒩, price_sell, palette = :tab10, label = "price sell")
        # push!(plots, p1)
        display(plot(p1, p2, layout =(2,1), title = "$n.id scenario $ω"))
    end
end


function plot_openloop(ocp::Model, nodes::Array, mpc::AbstractMPC)
    plot_energy_balance_openloop(ocp, nodes, mpc)

    for n in nodes
        plot_openloop(ocp, n, mpc)
    end
end

function plot_energy_balance(nodes, states, N)
    cumulative_output = zeros(N)
    cumulative_input = zeros(N)
    p = plot(title = "energybalance", ylimits=(-15, 15))
    for ((name,), node) in enamerate(nodes)
        if typeof(node) <: GasTurbine
            gt_power = states[name]["power"]
            plot!(1:N, cumulative_output + gt_power, fillrange = cumulative_output, fillalpha = 0.2, palette = :tab10,label = name)
            cumulative_output += gt_power
        elseif typeof(node) <: WindTurbineSDE
            wind_power = states[name]["power"]
                plot!(1:N, cumulative_output + wind_power , fillrange = cumulative_output, fillalpha = 0.2, palette = :tab10,label =  name)
                cumulative_output += wind_power
        elseif typeof(node) <: Demand
            power = states[name]["P"]
            plot!(1:N, -cumulative_input - power, fillrange = -cumulative_input, fillalpha = 0.2, palette = :tab10,label =  name)
            cumulative_input += power

        elseif typeof(node) <: Battery
            power_in = states[name]["Pᵢₙ"]
            plot!(1:N, -cumulative_input - power_in, fillrange = -cumulative_input, fillalpha = 0.2,palette = :tab10,label =  "$name p_in")
            cumulative_input += power_in

            power_out = states[name]["Pₒᵤₜ"]
            plot!(1:N, cumulative_output + power_out, fillrange = cumulative_output, fillalpha = 0.2, palette = :tab10,label =  "$name p_out")
            cumulative_output += power_out
        elseif typeof(node) <: SimpleMarket
            power_sell = states[name]["power_sell"]
            power_buy = states[name]["power_buy"]

            plot!(1:N, -cumulative_input - power_sell, fillrange = -cumulative_input, fillalpha = 0.2,palette = :tab10,label =  "$name sell")
            cumulative_input += power_sell

            plot!(1:N, cumulative_output + power_buy, fillrange = cumulative_output, fillalpha = 0.2, palette = :tab10,label =  "$name buy")
            cumulative_output += power_buy
        end
    end

    return(p, cumulative_input, cumulative_output)
end


function plot_node(name, node::GasTurbine, states, inputs, N::Int)
    p1 = plot(1:N, inputs[name]["load"], palette = :tab10, label = "load")
    p1 = plot!(1:N, states[name]["η"], palette = :tab10, label = "efficiency", ylimits=(0,1))
    p2 = plot(1:N, states[name]["power"], palette = :tab10, label = "power", ylimits=(0,node.Pₙₒₘ))

    return([p1,p2])
    # display(plot!(p1, p2, layout =(2,1), title="closed loop $name"))
end


function plot_node(name, node::Battery, states, inputs, N::Int)
    #battery soc
    p1 = plot(1:N, states[name]["soc"], palette = :tab10, label = "SOC", ylimits=(0,1))
    # p2 = plot(1:N, inputs[name]["charge"], palette = :tab10, label = "charge", ylimits=(-1,1))
    # p2 = plot!(1:N, -inputs[name]["discharge"], palette = :tab10, label = "discharge")
    p3 = plot(1:N, states[name]["Pᵢₙ"], palette = :tab10, label = "power in", ylimits=(-node.rate,node.rate))
    p3 = plot!(1:N, -states[name]["Pₒᵤₜ"], palette = :tab10, label = "power out")
    
    return([p1, p3])
    # display(plot(p1, p2, p3, layout =(3,1), title="closed loop $name"))
end

function plot_node(name, node::WindTurbineSDE, states, inputs, N::Int)
    p1 = plot(1:N, states[name]["power"], fillrange = zero(N), fillalpha = 0.2, palette = :tab10, label = "wind power",  ylimits=(-node.Pₘₐₓ,node.Pₘₐₓ))
    p1 = plot!(1:N, -inputs[name]["curt"]*node.Pₘₐₓ, fillrange = zero(N), fillalpha = 0.2, palette = :tab10, label = "power curtailed")
    p1 = plot!(;title="closed loop $name")
    return([p1])
    # display(p1)
end

function plot_node(name, node::Demand, states, inputs, N::Int)
    p = plot(1:N, states[name]["P"], palette = :tab10, fillrange = zero(N), fillalpha = 0.2, label = "Power in", title="demand")
    return([p])
end

function plot_node(name, node::SimpleMarket, states, inputs, N::Int)
    p1 = plot(1:N, -states[name]["power_buy"], palette = :tab10, fillrange = zero(N), fillalpha = 0.2, label = "buy", ylimits=(-node.capacity,node.capacity))
    p1 = plot!(1:N, states[name]["power_sell"], palette = :tab10, fillrange = zero(N), fillalpha = 0.2, label = "sell", ylimits=(-node.capacity,node.capacity))
    return([p1])
end

function process_history(history::Dict, data::Dict, nodes::NamedVector, N_sim::Int, filename::String, save_results=true)
    # get all x trajectories from feedback
    states = Dict(k => Dict() for k in collect(Iterators.flatten(names(nodes))))
    inputs = Dict(k => Dict() for k in collect(Iterators.flatten(names(nodes))))

    for n in nodes.array
        x_idx, u_idx = define_variable_indeces!(n)
        for xid in x_idx
            states[n.id][xid] = []
        end
        for uid in u_idx
            inputs[n.id][uid] = []
        end
    end

    for k in range(1,N_sim)

        action = history["action"][k]
        for ((name,), node_action) in enamerate(action)
            if typeof(node_action) == GTAction
                push!(inputs[name]["load"], node_action.load) 
            elseif typeof(node_action) == BatteryAction
                push!(inputs[name]["charge"], node_action.charge)
                push!(inputs[name]["discharge"], node_action.discharge)
            elseif typeof(node_action) == WindAction
                push!(inputs[name]["curt"], node_action.curtailment)
            elseif typeof(node_action) == SimpleMarketAction
                push!(inputs[name]["sell"], node_action.sell)
                push!(inputs[name]["buy"], node_action.buy)
                # @error("action node not accounted for")
            end
        end

        feedback = history["feedback"][k]
        for ((name,), node_feedback) in enamerate(feedback)
            if typeof(node_feedback) == GTFeedback
                push!(states[name]["η"], node_feedback.η)
                push!(states[name]["power"], node_feedback.power)
                push!(states[name]["CO₂"], node_feedback.co2)
            elseif typeof(node_feedback) == BatteryFeedback
                push!(states[name]["soc"], node_feedback.soc)
                push!(states[name]["Pₒᵤₜ"], node_feedback.power_out)
                push!(states[name]["Pᵢₙ"], node_feedback.power_in)

            elseif typeof(node_feedback) == WindFeedback
                push!(states[name]["power"], node_feedback.power)
            elseif typeof(node_feedback) == DemandFeedback
                push!(states[name]["P"], node_feedback.power)
            elseif typeof(node_feedback) == SimpleMarketFeedback
                push!(states[name]["power_buy"], node_feedback.power_buy)
                push!(states[name]["power_sell"], node_feedback.power_sell)
            else
                @error("feedback node not accounted for")
            end
            
        end
    end

    # ToDo because of the market timing we need to add initial market conditions
    simplemarkets_actions = filter(x -> isa(x, SimpleMarketAction),  history["action"][1])
    if length(simplemarkets_actions) > 0
        for sma in simplemarkets_actions
            states[sma.id]["power_buy"] = vcat([0.0], states[sma.id]["power_buy"][1:end-1])
            states[sma.id]["power_sell"] = vcat([0.0], states[sma.id]["power_sell"][1:end-1])
        end
    end


    if save_results==true
        if isnothing(filename)
            #path = "./Results/"*data["start"]*"_"*data["stop"]*".jld"
            path = "./Results/"* Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")*".jld2"
            save(path,  "N_sim", N_sim, "nodes", nodes, "states", states, "inputs", inputs, "data", data)
        end
        d =  Hour(data["stop"]-data["start"]).value/24
        path = "./Results/"* Dates.format(data["start"], "yyyy-mm-ddTHH:MM:SS")*"/$d days/" *filename*".jld2"
            save(path, "N_sim", N_sim, "nodes", nodes, "states", states, "inputs", inputs, "data", data, "openloop", history["openloop"])
    end

    return states, inputs

end

function plot_closed_loop(nodes::NamedVector, history::Dict, data::Dict, filename::String)
    N_sim = length(history["demand"])
    states, inputs = process_history(history, data, nodes, N_sim, filename)

    # Call all the plotting functions here
    p, cumulative_input, cumulative_output = plot_energy_balance(nodes, states, N_sim)
    p = plot!(title="energy balance")
    display(p)
    p1 = plot(cumulative_input)
    p1 = plot!(cumulative_output, label="cumulative output")
    display(p1)

    for ((name,), node) in enamerate(nodes)
        p = plot_node(name, node, states, inputs, N_sim)
        display(plot(p..., layout =(length(p),1), title="closed loop $name"))
    end

    return states, inputs

end

function plot_history(path::String)
    d = load(path)
    states = d["states"]
    inputs = d["inputs"]
    nodes = d["nodes"]
    N_sim = d["N_sim"]
    # Call all the plotting functions here
    display(plot(d["data"]["power_pot"], title="windpower"))
    p, _,_ = plot_energy_balance(nodes, states, N_sim)
    display(p)

    for ((name,), node) in enamerate(nodes)
        p = plot_node(name, node, states, inputs, N_sim)
        display(plot(p..., layout =(length(p),1), title="closed loop $name"))
    end
end


function plot_history(path1::String, path2::String)
    d1 = load(path1)
    states1 = d1["states"]
    inputs1 = d1["inputs"]
    nodes1 = d1["nodes"]
    N_sim1 = d1["N_sim"]

    d2 = load(path2)
    states2 = d2["states"]
    inputs2 = d2["inputs"]
    nodes2 = d2["nodes"]
    N_sim2 = d2["N_sim"]

    # Call all the plotting functions here
    display(plot(d1["data"]["power_pot"], title="windpower"))

    # path1
    p, input, output = plot_energy_balance(nodes1, states1, N_sim1)
    plot_name1 = split(path1, "/")[end]
    plot_name1 = split(plot_name1, ".")[1]
    plot_name1 = split(plot_name1, "_")[1]
    p = plot!(title=plot_name1)
    display(p)

    imbalance = sum(max.(0, input.-output))
    @info("imbalance $plot_name1 : $imbalance")
    plot(input, label="cumulative_input")
    display(plot!(output,label="cumulative_output"))

    # path2
    p, input, output = plot_energy_balance(nodes2, states2, N_sim2)
    plot_name2 = split(path2, "/")[end]
    plot_name2 = split(plot_name2, ".")[1]
    plot_name2 = split(plot_name2, "_")[1]
    p = plot!(title=plot_name2)
    display(p)

    imbalance = sum(max.(0, input.-output))
    @info("imbalance $plot_name2 : $imbalance")
    plot(input, label="cumulative_input")
    display(plot!(output,label="cumulative_output"))


    for ((name,), node) in enamerate(nodes1)
        p = plot_node(name, node, states1, inputs1, N_sim1)
        display(plot(p..., layout =(length(p),1), title="$name $plot_name1"))

        p = plot_node(name, node, states2, inputs2, N_sim2)
        display(plot(p..., layout =(length(p),1), title="$name $plot_name2"))
    end
end