function cost_balance(ocp::Model, mpc::DeterministicMPC)
    V = 0.0
    𝒩 = mpc.𝒩
    for k in 𝒩
        l = 0.0
        l += mpc.weights["balance"] * ocp[:balance][k]
        V += l
    end
    return V
end

function cost_balance(ocp::Model, mpc::RiskMeasureMPC, ω::Int, k::Int)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω
    lω = mpc.weights["balance"] * ocp[:balance][ω, k]
    return lω
end

function get_cost(ocp::Model, node::GasTurbine, mpc::DeterministicMPC)
    V = 0.0
    𝒩 = mpc.𝒩
    for k in 𝒩
        l = 0.0
        # l += mpc.weights["Δgt"] * (ocp[:u][node.id, "load", k+1]-ocp[:u][node.id, "load", k])^2
        l += mpc.weights["Δgt"] * ocp[:u][node.id, "load", k]^2
        V += l
    end
    return V
end

function get_cost(ocp::Model, node::GasTurbine, mpc::RiskMeasureMPC, ω::Int, k::Int)
    V = 0.0
    𝒩 = mpc.𝒩

    # l += mpc.weights["Δgt"] * (ocp[:u][node.id, "load", k+1]-ocp[:u][node.id, "load", k])^2
    lω = mpc.weights["Δgt"] * ocp[:u][node.id, "load", ω, k]^2
    return lω
end



function get_cost(ocp::Model, node::Battery, mpc::DeterministicMPC)
    V = 0.0
    𝒩 = mpc.𝒩

    for k in 𝒩
        l = 0.0
        hubber_in = sqrt(0.05 +(ocp[:x][node.id, "Pᵢₙ", k]/node.rate)^2) - 0.05
        hubber_out = sqrt(0.05 +(ocp[:x][node.id, "Pₒᵤₜ", k]/node.rate)^2) - 0.05
        l += mpc.weights["bat"] * (hubber_in + hubber_out)

        #penalty to prevent charging and discharging at the same time
        l += (mpc.weights["bat_chargedischarge"] *  (ocp[:x][node.id, "Pᵢₙ", k] 
        * ocp[:x][node.id, "Pₒᵤₜ", k]))

        # Battery SOC penalty 10-90
        l += (mpc.weights["bat_min"] * (ocp[:σ][node.id, "socₘᵢₙ", k]))
        l += (mpc.weights["bat_max"] * (ocp[:σ][node.id, "socₘₐₓ", k]))

        V += l
    end

    #Battery terminal cost constraint
    l = 0.0 
    l = mpc.weights["bat_terminal"] *(ocp[:x][node.id, "soc", k] - 0.5)^2
    l += mpc.weights["bat_terminal"] *(ocp[:x][node.id, "soc", k] - 0.5)^2
    
    return V
end

function get_cost(ocp::Model, node::Battery, mpc::RiskMeasureMPC, ω::Int, k::Int)
    V = 0.0
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω

    lω = 0.0
    hubber_inω = sqrt(0.05 +(ocp[:x][node.id, "Pᵢₙ", ω, k]/node.rate)^2) - 0.05
    hubber_outω = sqrt(0.05 +(ocp[:x][node.id, "Pₒᵤₜ", ω, k]/node.rate)^2) - 0.05
    lω += mpc.weights["bat"] * (hubber_inω + hubber_outω)

    #penalty to prevent charging and discharging at the same time
    lω += (mpc.weights["bat_chargedischarge"] *  (ocp[:x][node.id, "Pᵢₙ", ω, k] 
    * ocp[:x][node.id, "Pₒᵤₜ", ω, k]))

    # Battery SOC penalty 10-90
    lω += (mpc.weights["bat_min"] * (ocp[:σ][node.id, "socₘᵢₙ", ω, k]))
    lω += (mpc.weights["bat_max"] * (ocp[:σ][node.id, "socₘₐₓ", ω, k]))

    if k == last(𝒩)
        lω += mpc.weights["bat_terminal"] *(ocp[:x][node.id, "soc", ω, k] - 0.5)^2
    end
    
    return lω
end

function get_cost(ocp::Model, node::WindTurbineSDE, mpc::DeterministicMPC)
    V = 0.0
    𝒩 = mpc.𝒩
    
    return V
end

function get_cost(ocp::Model, node::WindTurbineSDE, mpc::RiskMeasureMPC, ω::Int, k::Int)
    lω = 0.0
    𝒩 = mpc.𝒩
    return lω
end


function get_cost(ocp::Model, node::Demand, mpc::DeterministicMPC)
    lω = 0.0
    𝒩 = mpc.𝒩
    return lω
end

function get_cost(ocp::Model, node::Demand, mpc::RiskMeasureMPC, ω::Int, k::Int)
    lω = 0.0
    𝒩 = mpc.𝒩
    return lω
end

function get_cost(ocp::Model, node::SimpleMarket, mpc::DeterministicMPC)
    lω = 0.0
    𝒩 = mpc.𝒩
    for k in 𝒩
        lω += mpc.weights["market"] * (ocp[:price_buy][node.id, k]) * ocp[:x][node.id, "power_buy", k]
        lω -= mpc.weights["market"] * (ocp[:price_sell][node.id, k]) * ocp[:x][node.id, "power_sell", k]
    end
    return lω
end

function get_cost(ocp::Model, node::SimpleMarket, mpc::RiskMeasureMPC, ω::Int, k::Int)
    lω = 0.0
    lω += mpc.weights["market"] * (ocp[:price_buy][node.id, k]) * ocp[:x][node.id, "power_buy", ω, k]
    lω -= mpc.weights["market"] * (ocp[:price_sell][node.id, k]) * ocp[:x][node.id, "power_sell", ω, k]
    return lω
end


# function cost(OCP::Model, node::Demand, k::Int, ω::Int, mpc::Scenario_SMPC)
#     l = 0
#     l = (mpc.weights["loadshift"] * OCP[:u][k, node.u[1], ω]^2)
#     l += (mpc.weights["demand_ΔE"] * OCP[:x][k, node.x[2], ω])

#     return l
# end

