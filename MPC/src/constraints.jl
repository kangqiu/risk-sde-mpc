function model_constraints!(ocp::Model, n::GasTurbine, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    

    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "η", k] == 
    (ocp[:u][n.id, "load", k]*n.ηₛ)
    / (ocp[:u][n.id, "load", k] * n.ηₛ + (1-n.ηₛ)*((1-n.Fᵪ)+n.Fᵪ*ocp[:u][n.id, "load", k]^2))
    )

    #GT power output
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "power", k] == n.Pₙₒₘ*(ocp[:x][n.id, "η", k]/n.ηₛ ))

    #GT CO2 output
    #CH₄ + 2O₂ → 1CO₂ + 2H₂O 
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "CO₂", k] == ocp[:x][n.id, "power", k] * mpc.Δt /n.LHV * (n.M_CO₂/n.M_CH₄)
    )

end

function model_constraints!(ocp::Model, n::GasTurbine, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "η", ω, k] == 
    (ocp[:u][n.id, "load", ω,  k] * n.ηₛ)
    / (ocp[:u][n.id, "load", ω, k] * n.ηₛ + (1-n.ηₛ)*((1-n.Fᵪ)+n.Fᵪ*ocp[:u][n.id, "load", ω, k]^2))
    )

    #GT power output
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "power", ω, k] == n.Pₙₒₘ*(ocp[:x][n.id, "η", ω, k]/n.ηₛ ))

    #GT CO2 output
    #CH₄ + 2O₂ → 1CO₂ + 2H₂O 
    @constraint(ocp, [k ∈ 𝒩,  ω ∈ Ω],
    ocp[:x][n.id, "CO₂", ω, k] == ocp[:x][n.id, "power", ω, k] * mpc.Δt /n.LHV * (n.M_CO₂/n.M_CH₄)
    )

end

function model_constraints!(ocp::Model, n::Battery, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    # Battery SOC
    @constraint(ocp, [k ∈ 𝒩[1:end-1]],
    ocp[:x][n.id, "soc", k+1] == ocp[:x][n.id, "soc", k]+1/n.capacity * mpc.Δt * n.rate * (n.η*ocp[:u][n.id, "charge", k]-n.η^(-1)*
    ocp[:u][n.id, "discharge", k]) 
    )

    # Battery power constraints
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "Pₒᵤₜ", k] == n.rate * n.η^(-1)*ocp[:u][n.id, "discharge", k]) 
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "Pᵢₙ", k] == n.rate*n.η*ocp[:u][n.id, "charge", k]) 

end

function model_constraints!(ocp::Model, n::Battery, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω
    # Battery SOC
    @constraint(ocp, [k ∈ 𝒩[1:end-1], ω ∈ Ω],
    ocp[:x][n.id, "soc", ω, k+1] == ocp[:x][n.id, "soc", ω, k]+1/n.capacity * mpc.Δt * n.rate * (n.η*ocp[:u][n.id, "charge", ω, k]-n.η^(-1)*
    ocp[:u][n.id, "discharge", ω, k]) 
    )

    # Battery power constraints
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "Pₒᵤₜ", ω, k] == n.rate * n.η^(-1)*ocp[:u][n.id, "discharge", ω, k]) 
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "Pᵢₙ", ω, k] == n.rate*n.η*ocp[:u][n.id, "charge", ω, k]) 

end




function model_constraints!(ocp::Model, n::WindTurbineSDE, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "power", k] == n.Pₘₐₓ * (ocp[:windpower][n.id, k] - ocp[:u][n.id, "curt", k]) )
    
end


function model_constraints!(ocp::Model, n::WindTurbineSDE, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "power", ω, k] == n.Pₘₐₓ * (ocp[:windpower][n.id, ω, k] - ocp[:u][n.id, "curt", ω, k]) )
    
end


function model_constraints!(ocp::Model, n::Demand, mpc::AbstractMPC)
    𝒩 = mpc.𝒩
end


function model_constraints!(ocp::Model, n::SimpleMarket, mpc::AbstractMPC)
    𝒩 = mpc.𝒩
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "power_buy", k] == n.capacity * ocp[:u][n.id, "buy", k]) 
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "power_sell", k] == n.capacity * ocp[:u][n.id, "sell", k]) 
end

function model_constraints!(ocp::Model, n::SimpleMarket, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω = 1:mpc.model.ω
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "power_buy", ω, k] == n.capacity * ocp[:u][n.id, "buy", ω, k]) 
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "power_sell", ω, k] == n.capacity * ocp[:u][n.id, "sell", ω, k]) 
end



function state_constraints!(ocp::Model, n::GasTurbine, mpc::AbstractMPC)
    𝒩 = mpc.𝒩

end




function state_constraints!(ocp::Model, n::Battery, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    #SOC constraint
    @constraint(ocp, [k ∈ 𝒩], 0.0<= ocp[:x][n.id, "soc", k] <= 1.0)

    #SOC initial value constraint
    @constraint(ocp, 
    ocp[:x][n.id, "soc", 1] == ocp[:x₀][n.id, "soc"]
    )

    # slack variables
    @constraint(ocp, [k ∈ 𝒩], ocp[:σ][n.id, "socₘᵢₙ", k] >= 0.0)
    @constraint(ocp, [k ∈ 𝒩], ocp[:σ][n.id, "socₘᵢₙ", k] >= (n.socₘᵢₙ - ocp[:x][n.id, "soc", k])*10) #normalizes the slack variable between 0 and 1
    @constraint(ocp, [k ∈ 𝒩], ocp[:σ][n.id, "socₘₐₓ", k] >= 0.0)
    @constraint(ocp, [k ∈ 𝒩], ocp[:σ][n.id, "socₘₐₓ", k] >= (ocp[:x][n.id, "soc", k] - n.socₘₐₓ)*10)

end

function state_constraints!(ocp::Model, n::Battery, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    #SOC constraint
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], 0.0 <= ocp[:x][n.id, "soc", ω, k] <= 1.0)

    # #SOC initial value constraint
    @constraint(ocp, [ω ∈ Ω], 
    ocp[:x][n.id, "soc", ω, 1] == ocp[:x₀][n.id, "soc"]
    )

    # slack variables
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:σ][n.id, "socₘᵢₙ", ω, k] >= 0.0)
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:σ][n.id, "socₘᵢₙ", ω, k] >= (n.socₘᵢₙ - ocp[:x][n.id, "soc", ω, k])*10) #normalizes the slack variable between 0 and 1
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:σ][n.id, "socₘₐₓ", ω, k] >= 0.0)
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:σ][n.id, "socₘₐₓ", ω, k] >= (ocp[:x][n.id, "soc", ω, k] - n.socₘₐₓ)*10)

end


function state_constraints!(ocp::Model, n::WindTurbineSDE, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    @constraint(ocp, [k ∈ 𝒩],
    0.0 <= ocp[:x][n.id, "power", k] <= n.Pₘₐₓ
    )
    
end


function state_constraints!(ocp::Model, n::WindTurbineSDE, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    0.0 <= ocp[:x][n.id, "power", ω, k] <= n.Pₘₐₓ
    )
    
end

function state_constraints!(ocp::Model, n::Demand, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    @constraint(ocp, [k ∈ 𝒩],
    ocp[:x][n.id, "P", k] == ocp[:demand][n.id, k]
    )

end

function state_constraints!(ocp::Model, n::Demand, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
    ocp[:x][n.id, "P", ω, k] == ocp[:demand][n.id, k]
    )

end

function state_constraints!(ocp::Model, n::SimpleMarket, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    # initial market constraints
    @constraint(ocp,
    ocp[:x][n.id, "power_buy", 1] == ocp[:x₀][n.id, "power_buy"]
    )
    @constraint(ocp,
    ocp[:x][n.id, "power_sell", 1] == ocp[:x₀][n.id, "power_sell"]
    )

    # market capacity constraints
    @constraint(ocp, [k ∈ 𝒩],
        ocp[:x][n.id, "power_buy", k] <= n.capacity 
    )
    @constraint(ocp, [k ∈ 𝒩],
        ocp[:x][n.id, "power_buy", k] >= 0.0 
    )
    @constraint(ocp, [k ∈ 𝒩],
        ocp[:x][n.id, "power_sell", k] <= n.capacity 
    )
    @constraint(ocp, [k ∈ 𝒩],
        ocp[:x][n.id, "power_sell", k] >= 0.0 
    )
end

function state_constraints!(ocp::Model, n::SimpleMarket, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω
    # initial market constraints
    @constraint(ocp, [ω ∈ Ω],
    ocp[:x][n.id, "power_buy", ω, 1] == ocp[:x₀][n.id, "power_buy"]
    )
    @constraint(ocp, [ω ∈ Ω],
    ocp[:x][n.id, "power_sell", ω, 1] == ocp[:x₀][n.id, "power_sell"]
    )

    # market capacity constraints
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
        ocp[:x][n.id, "power_buy", ω, k] <= n.capacity 
    )
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
        ocp[:x][n.id, "power_buy", ω, k] >= 0.0 
    )
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
        ocp[:x][n.id, "power_sell", ω, k] <= n.capacity 
    )
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω],
        ocp[:x][n.id, "power_sell", ω, k] >= 0.0 
    )
end

function input_constraints!(ocp::Model, n::GasTurbine, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    #load constraint
    @constraint(ocp, [k ∈ 𝒩], n.loadₘᵢₙ <= ocp[:u][n.id, "load", k] <= 1.0)

end

function input_constraints!(ocp::Model, n::GasTurbine, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    #load constraint
    @constraint(ocp, [k ∈ 𝒩,  ω ∈ Ω], n.loadₘᵢₙ <= ocp[:u][n.id, "load", ω, k] <= 1.0)

    for ω ∈ Ω[1:end-1]
        @constraint(ocp, ocp[:u][n.id, "load", ω, 1] == ocp[:u][n.id, "load", ω+1, 1])
    end

end

function input_constraints!(ocp::Model, n::Battery, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    #charge discharge constraints
    @constraint(ocp, [k ∈ 𝒩], 0.0 <= ocp[:u][n.id, "charge", k] <= 1.0)
    @constraint(ocp, [k ∈ 𝒩], 0.0 <= ocp[:u][n.id, "discharge", k] <= 1.0)

end

function input_constraints!(ocp::Model, n::Battery, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    # #charge discharge constraints
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], 0.0 <= ocp[:u][n.id, "charge", ω, k] <= 1.0)
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], 0.0 <= ocp[:u][n.id, "discharge", ω, k] <= 1.0)

        for ω ∈ Ω[1:end-1]
            @constraint(ocp, ocp[:u][n.id, "charge", ω, 1] == ocp[:u][n.id, "charge", ω+1, 1])
            @constraint(ocp, ocp[:u][n.id, "discharge", ω, 1] == ocp[:u][n.id, "discharge", ω+1, 1])
        end

end

function input_constraints!(ocp::Model, n::WindTurbineSDE, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    @constraint(ocp, [k ∈ 𝒩], 0.0 <= ocp[:u][n.id, "curt", k])
    @constraint(ocp, [k ∈ 𝒩], ocp[:u][n.id, "curt", k] <= 1) # ocp[:windpower][n.id, k])

end

function input_constraints!(ocp::Model, n::WindTurbineSDE, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], 0.0 <= ocp[:u][n.id, "curt", ω, k])
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:u][n.id, "curt",ω, k] <= 1)


    # for k ∈ 𝒩
        for ω ∈ Ω[1:end-1]
            @constraint(ocp, ocp[:u][n.id, "curt", ω, 1] == ocp[:u][n.id, "curt", ω+1, 1])
        end
    # end

end

function input_constraints!(ocp::Model, n::Demand, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    @constraint(ocp, [k ∈ 𝒩], ocp[:u][n.id, "ΔP", k] ==0.0)
end


function input_constraints!(ocp::Model, n::Demand, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:u][n.id, "ΔP", ω, k] == 0.0)
end

function input_constraints!(ocp::Model, n::SimpleMarket, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩
    @constraint(ocp, [k ∈ 𝒩], ocp[:u][n.id, "buy", k] <= 1.0)
    @constraint(ocp, [k ∈ 𝒩], ocp[:u][n.id, "buy", k] >= 0.0)
    @constraint(ocp, [k ∈ 𝒩], ocp[:u][n.id, "sell", k] <= 1.0)
    @constraint(ocp, [k ∈ 𝒩], ocp[:u][n.id, "sell", k] >= 0.0)
end

function input_constraints!(ocp::Model, n::SimpleMarket, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:u][n.id, "buy", ω, k] <= 1.0)
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:u][n.id, "buy", ω, k] >= 0.0)
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:u][n.id, "sell", ω, k] <= 1.0)
    @constraint(ocp, [k ∈ 𝒩, ω ∈ Ω], ocp[:u][n.id, "sell", ω, k] >= 0.0)

    #non anticipativity
    for ω ∈ Ω[1:end-1]
        @constraint(ocp, ocp[:u][n.id, "buy", ω, 2] == ocp[:u][n.id, "buy", ω+1, 2])
        @constraint(ocp, ocp[:u][n.id, "sell", ω, 2] == ocp[:u][n.id, "sell", ω+1, 2])
    end
end



function add_energybalance!(ocp::Model, nodes::Vector, mpc::DeterministicMPC)
    𝒩 = mpc.𝒩

    σ = 0
    for n in nodes
        if typeof(n) <: WindTurbineSDE
                σ = σ .+ ocp[:x][n.id, "power", 𝒩]
        elseif typeof(n) <: GasTurbine
            σ = σ .+ ocp[:x][n.id, "power", 𝒩]
        elseif typeof(n) <: Battery
            σ = σ .+ ocp[:x][n.id, "Pₒᵤₜ", 𝒩]
            σ = σ .- ocp[:x][n.id, "Pᵢₙ", 𝒩]
        elseif typeof(n) <: Demand
            σ = σ .- ocp[:x][n.id, "P", 𝒩]
        elseif typeof(n) <: SimpleMarket
            σ = σ .- ocp[:x][n.id, "power_sell", 𝒩]
            σ = σ .+ ocp[:x][n.id, "power_buy", 𝒩]
        end
    end

    @constraint(ocp, [k ∈ 𝒩], ocp[:balance][k] >= σ[k])
    @constraint(ocp, [k ∈ 𝒩], -ocp[:balance][k] <= σ[k])

    # @constraint(ocp, [k ∈ 𝒩], σ[k] == 0.0)
    # @constraint(ocp, [k ∈ 𝒩], ocp[:balance][k] == 0.0)
end

function add_energybalance!(ocp::Model, nodes::Vector, mpc::RiskMeasureMPC)
    𝒩 = mpc.𝒩
    Ω =  1:mpc.model.ω

    for ω ∈ Ω
        σ = 0.0
        for n in nodes
            if typeof(n) <: WindTurbineSDE
                    σ = σ .+ ocp[:x][n.id, "power", ω, 𝒩]
            elseif typeof(n) <: GasTurbine
                σ = σ .+ ocp[:x][n.id, "power", ω, 𝒩]
            elseif typeof(n) <: Battery
                σ = σ .+ ocp[:x][n.id, "Pₒᵤₜ", ω, 𝒩]
                σ = σ .- ocp[:x][n.id, "Pᵢₙ", ω, 𝒩]
            elseif typeof(n) <: Demand
                σ = σ .- ocp[:x][n.id, "P", ω, 𝒩]
            elseif typeof(n) <: SimpleMarket
                σ = σ .- ocp[:x][n.id, "power_sell", ω, 𝒩]
                σ = σ .+ ocp[:x][n.id, "power_buy", ω, 𝒩]
            end
        end

        # equality balance
        @constraint(ocp, [k ∈ 𝒩], ocp[:balance][ω, k] >= σ[k])
        @constraint(ocp, [k ∈ 𝒩], -ocp[:balance][ω, k] <= σ[k])

        #inequality balance, meaning we can tolerate more production but not more demand
        # @constraint(ocp, [k ∈ 𝒩], ocp[:balance][ω, k] >= -σ[k])
        # @constraint(ocp, [k ∈ 𝒩], ocp[:balance][ω, k] >= 0)
        # @constraint(ocp, [k ∈ 𝒩], σ[k] >= 0)

    end
end