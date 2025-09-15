"""
    GT_partload_efficiency(system_data::Dict, gt_load::Float64)

Return the calculated partload efficiency ηₚₗ of a single cycle gas turbine. 

Based on the GE LM2500 aeroderivative gas turbines. Equations from [1] F. Haglind, B. Elmegaard,
Methodologies for predicting the part-load performance of aero-derivative gas turbines,
Energy, 2009,https://doi.org/10.1016/j.energy.2009.06.042.

# Arguments
- `u::Float64`: gas turbine load setpoint. u ∈ [0.35, 1.00] 
- `ηₛ::Float64=0.38`: nominal single cycle efficiency.
- `Fᵪ::Float64=0.43`: copper losses [1].
- `LHV::Float64=55.5`: lower heating value of natural gas.
"""
function efficiency(u::Real, n::GasTurbine)
    ηₚₗ = (u*n.ηₛ)/(u*n.ηₛ + (1-n.ηₛ)*((1-n.Fᵪ)+n.Fᵪ*u^2))
    return ηₚₗ
end

"""
    GT_power_output(system_data::Dict, gt_load::Float64)

Return the power output P of single cycle gas turbine.

Based on the GE LM2500 aeroderivative gas turbines. Assumed quasistatic as GTs have a quick settling time ~60s.

# Arguments
- `u::Float64`: gas turbine load setpoint. u ∈ [0.35, 1.00]
- `η::Float64`: gas turbine partload efficiency
- `Pₙₒₘ::Float64=34.0`: nominal capacity in MW.
- `ηₛ::Float64=0.34`: nominal single cycle efficiency.
"""

function GT_power(η::Real, n::GasTurbine)
    P = n.Pₙₒₘ*(η/n.ηₛ)
    return P
end

function co2_emissions(power::Real, n::GasTurbine, Δt::Real)
    return co₂ = power * Δt /n.LHV * (n.M_CO₂/n.M_CH₄)
end

"""
    battery(x::Float64)
Return the battery state of charge
"""

function stateofcharge(power_in::Real, power_out::Real, x_soc::Real, Δt::Real, n::Battery)

    soc = x_soc+ 1/n.capacity * Δt * ( n.η * power_in -n.η^(-1)*
    power_out) 

    soc = clamp(soc, 0,1)

    return soc
end

function battery_power_out(discharge::Real, n::Battery)
    P = n.rate * n.η ^-1 *discharge
    return P
end

function battery_power_in(charge::Real, n::Battery)
    P = n.rate*n.η*charge
    return P
end

function wind_curve(windspeed)
        windcurve = 1 -1 /(1 + ((windspeed +17)/24)^20) 
        return windcurve
end

function simple_market(buy::Real, sell::Real, n::SimpleMarket)
    power_bought = clamp(buy * n.capacity, 0, n.capacity)
    power_sold = clamp(sell * n.capacity, 0, n.capacity)
    return power_bought, power_sold
end

function simulate_system!(data::Dict, nodes::NamedArray, actions::NamedArray, feedback::NamedArray, k::Int,
    mpc::AbstractMPC)

    # Update feedback and actions based on the current state of the system
    cumulative_input = 0.0
    cumulative_output = 0.0
    wind_output = 0.0
    wind_output_pot = 0.0
    for node in nodes.array
        if node isa GasTurbine
            η = efficiency(actions[node.id].load, node)
            power = GT_power(η, node)
            cumulative_output += power
            feedback[node.id].η = η
            feedback[node.id].power = power
            feedback[node.id].co2 = co2_emissions(power, node, mpc.Δt)
        elseif node isa Battery
            power_out = battery_power_out(actions[node.id].discharge, node)
            power_in = battery_power_in(actions[node.id].charge, node)
            feedback[node.id].soc = stateofcharge(power_in, power_out,
                                                      feedback[node.id].soc, mpc.Δt, node)
            feedback[node.id].power_out = power_out
            feedback[node.id].power_in = power_in

            cumulative_output += power_out
            cumulative_input += power_in
        elseif node isa WindTurbineSDE
            windpower = data["power_pot"][k]
            power = node.Pₘₐₓ * (windpower - actions[node.id].curtailment)
            feedback[node.id].power = clamp(power, 0.0, node.Pₘₐₓ)
            feedback[node.id].power_pot = windpower
            feedback[node.id].v = data["v_wind"][k]
            feedback[node.id].r = data["r_wind"][k]
            feedback[node.id].q = data["q_wind"][k]

            wind_output += power
            wind_output_pot += windpower
        elseif node isa Demand
            feedback[node.id].power = data["demand"][k]
            cumulative_input += data["demand"][k]

        elseif node isa SimpleMarket
            feedback[node.id].power_buy, feedback[node.id].power_sell = simple_market(
                actions[node.id].buy, actions[node.id].sell, node)
            cumulative_input += feedback[node.id].power_sell
            cumulative_output += feedback[node.id].power_buy
        end
    end
    return feedback
end