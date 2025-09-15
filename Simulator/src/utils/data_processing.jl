"""
    get_historical_forecast_data(start::DateTime, stop::DateTime, mpc_horizon::Int; data_path::String="...")

Load and process historical forecast (NWP) data. Function is intended to be used before each simulation run. 
The forecasts are then provided to both simulation and MPC.

# Arguments
- `start::DateTime`: Simulation start time.
- `stop::DateTime`: Simulation stop time.
- `mpc_horizon::Int`: Model predictive control horizon in hours.
- `data_path::String`: Path to the JLD2 file containing forecast data.

# Returns
- `df_nwp`: DataFrame with NWP data and computed dnwp column with 10min time steps.
"""
function get_historical_forecast_data(start::DateTime, stop::DateTime, mpc_horizon::Int; data_path::String="Simulator/data/data_forecast.jld2")
    data = jldopen(data_path, "r")  
    data = data["forecast"]
    
    stop_datetime = stop + Hour(mpc_horizon)

    # get indices for data array (in order to interpolate later we add an hour before and after the simulation period)
    first_index = findlast(x -> x <= start - Hour(1), data[!, "time"])
    last_index = findfirst(x -> x >= stop_datetime + Hour(1), data[!, "time"])

    # retrieve nwp data
    df_nwp = data[first_index:last_index, :]

    # forecast update times for the simulation period
    forecast_update_times = collect(DateTime(2020,1,1, 0, 0, 0):Hour(6):DateTime(2025,01,01, 0, 0,0))

    index_start = findlast(x -> x <= df_nwp[1, "time"], forecast_update_times)
    index_stop = findfirst(x -> x >= df_nwp[end, "time"], forecast_update_times)
    forecast_update_times = forecast_update_times[index_start-1:index_stop]

    # calculate the dnwp, difference in weather prediction between updates every 6h
    dnwp = []
    index_nwp = findfirst(x-> x == forecast_update_times[1], data[!, "time"])
    nwp_0 = data[index_nwp, "forecast"]

    for ts in range(2,length(forecast_update_times))
        index_time = findfirst(x-> x == forecast_update_times[ts], data[!, "time"])
        push!(dnwp, data[index_time, "forecast"] - nwp_0)
        nwp_0 = data[index_time, "forecast"] 
    end
    
    #create interpolation object to downsample dnwp to times required for nwp DataFrame
    dnwp_x = (forecast_update_times[2:end] .- forecast_update_times[1])./Minute(10)
    itp = Interpolations.interpolate((dnwp_x,), dnwp, Gridded(Constant(Previous)))
    df_nwp_times = (df_nwp[!, "time"].-forecast_update_times[1])./Minute(10)

    dnwp_10min = itp.(df_nwp_times)

    df_nwp[!, "dnwp"] = dnwp_10min
    return df_nwp
end


"""
    get_interp_object_best_nwp(
        data::DataFrame, 
        start::DateTime, 
        stop::DateTime, 
        tspan::Tuple{<:Real, <:Real}
    )

Creates interpolation objects for wind forecast and forecast update increments (dnwp) over the simulation period. 
    This method uses the best available forecast.

# Arguments
- `data::DataFrame`: DataFrame containing NWP data with "forecast" and "dnwp" columns at 10-minute intervals.
- `start::DateTime`: Simulation start time.
- `stop::DateTime`: Simulation stop time.

# Returns
- `itp_nwp`: Interpolation object for wind forecast values. Input: time in hours since start of simulation.
- `itp_dnwp`: Interpolation object for dnwp values. Input: time in hours since start of simulation
"""

function get_interp_object_best_nwp(data::DataFrame, start::DateTime, stop::DateTime, Δt::Int,  N_mpc::Int)
    # last_index = findfirst(x -> x >= stop+Hour(1), data[!, "time"])
    times = data[:, "time"]
    forecasts = data[:, "forecast"]
    dforecasts = data[:, "dnwp"]

    # time grid in hours for interpolation (starts one hour before and ends one hour after the simulation period)
    time_grid_10 = -1:10/60:Hour((stop+Minute(Δt)*N_mpc + Hour(1))-start).value

    # fill nwp with best available forecast
    nwp = Float64[]
    dnwp = Float64[]
    for t in range(1, length(times))
        push!(nwp, forecasts[t][1])
        push!(dnwp, dforecasts[t][1])
    end

    itp_nwp = DataInterpolations.ConstantInterpolation(nwp, time_grid_10)
    itp_dnwp = DataInterpolations.ConstantInterpolation(dnwp, time_grid_10)

    return itp_nwp, itp_dnwp
end

"""
    get_mpc_forecast_interp(
        data::DataFrame,
        start::DateTime,
        stop::DateTime,
        T_mpc::Int,
        Δt::Int,
        k::Int
    )

Creates interpolation objects for wind forecast and forecast update increments (dnwp) for use in Model Predictive Control (MPC) over the MPC horizon.

# Arguments
- `data::DataFrame`: DataFrame containing NWP data with "forecast" and "dnwp" columns.
- `start::DateTime`: Simulation start time.
- `stop::DateTime`: Simulation stop time.
- `T_mpc::Int`: MPC horizon in hours.
- `Δt::Int`: Time step in minutes for the output time grid.
- `k::Int`: MPC step index.
"""

function get_mpc_forecast_interp(data::DataFrame, start::DateTime, stop::DateTime, T_mpc::Real,  Δt::Int, k::Int )

    start_mpc = start + Minute(Δt * (k)) # start time for the MPC step
    stop_mpc = start_mpc + Hour(T_mpc ) # stop time for the
    last_index = findfirst(x -> x >= stop_mpc + Hour(1), data[!, "time"])
    first_index = findlast(x -> x <= start_mpc - Hour(1), data[!, "time"])
    times = data[first_index:last_index, "time"]
    forecasts = data[first_index:last_index, "forecast"]
    dforecasts = data[first_index:last_index, "dnwp"]

    # time grid in hours for interpolation (starts one hour before and ends one hour after the simulation period)
    time_grid = -1:10/60:T_mpc+1

    # fill nwp with best available forecast
    nwp = Float64[]
    dnwp = Float64[]
    for t in range(1, length(times))
        push!(nwp, forecasts[t][1])
        push!(dnwp, dforecasts[t][1])
    end

    itp_nwp = DataInterpolations.ConstantInterpolation(nwp, time_grid)
    itp_dnwp = DataInterpolations.ConstantInterpolation(dnwp, time_grid)

    return itp_nwp, itp_dnwp
end

"""
    sample_Gaussian(mean::Array{<:Real}, std::Array{<:Real})

Generates samples from Gaussian distributions for each pair of mean and standard deviation using percentiles.

# Arguments
- `mean::Array{<:Real}`: Array of mean values for the Gaussian distributions from SDE in Lamperti domain.
- `std::Array{<:Real}`: Array of standard deviations for the Gaussian distributions from SDE in Lamperti domain.

# Returns
- `samples::Array`: Array where each element contains quantile samples from the corresponding Gaussian distribution in natural domain.

# Description
For each (mean, std) pair, this function creates a Normal distribution and computes quantiles at percentiles from 0.01 to 0.99, returning the samples as an array of arrays.
"""
function sample_power_Gaussian(mean_power::Array{<:Real}, std_power::Array{<:Real}, mean_wind::Array{<:Real}, std_wind::Array{<:Real},
     data::Dict, k::Int)
    percentiles = collect(0.01:0.01:0.99)
    # percentiles = [0.05, 0.5, 0.95]
    samples = []
    samples_wind = []
    for (i, (μ_p, σ_p, μ_v, σ_v)) in enumerate(zip(mean_power, std_power, mean_wind, std_wind))
        dist = Normal(μ_p, σ_p)
        dist_wind = Normal(μ_v, σ_v)
        s = quantile.(dist, percentiles)
        s_wind = quantile.(dist_wind, percentiles)

        s_n =[]
        s_wind_n = []
        for (j,α) in enumerate(percentiles)
            power_sample = lamperti_to_natural_transform_power(s[j], data["power_pot_l"][k+i], data["power_pot"][k+i])
            wind_sample = lamperti_to_natural_transform_wind(s_wind[j])

            #implement cutout
            # if wind_sample >= 25
            #     power_sample = 0
            # end

            push!(s_n, power_sample)
            push!(s_wind_n, wind_sample)

        end
        push!(samples, s_n)
        push!(samples_wind, s_wind_n)
    end
    return(samples, samples_wind)
end


"""
    get_samples(mpc::MPC.RiskMeasureMPC, df_nwp::DataFrame, data::Dict, k::Int)

Generates scenario samples for wind power output over the MPC horizon using NWP data and SDE simulation results.
This code works only for the non-graph based scenario representation and does not consider transition probabilities!
Alternative to scenario tree!

# Arguments
- `mpc::MPC.RiskMeasureMPC`: The risk-aware MPC object containing simulation parameters.
- `df_nwp::DataFrame`: DataFrame containing NWP data with "forecast" and "dnwp" columns.
- `data::Dict`: Dictionary with simulation and MPC configuration (e.g., start/stop times, horizon).
- `k::Int`: MPC step index.

# Returns
- `samples::Array`: Array of scenario samples for wind power output in the natural domain.

# Description
This function creates interpolation objects for wind forecast and its increments for the current MPC step, simulates the SDE to obtain mean and standard deviation for wind and power in the Lamperti domain, and generates quantile-based samples in the natural domain using Gaussian distributions.
"""

function get_samples(mpc::RiskMeasureMPC, df_nwp::DataFrame, data::Dict, k::Int)
    itp_mpc_nwp, itp_mpc_dnwp = get_mpc_forecast_interp(
        df_nwp, 
        data["start"], 
        data["stop"],
        data["T_mpc"],
        mpc.Δt,
        k);

    μ_power_L, σ_power_L, μ_wind_L, σ_wind_L = simulate_ode(
        data["T_mpc"], data, itp_mpc_nwp, itp_mpc_dnwp, mpc.Δt, k);

    
    samples, samples_wind = sample_power_Gaussian(μ_power_L, σ_power_L, μ_wind_L, σ_wind_L, data, k);

    return(samples, samples_wind)
    
end

function get_samples(mpc::DeterministicMPC, df_nwp::DataFrame, data::Dict, k::Int)
    itp_mpc_nwp, itp_mpc_dnwp = get_mpc_forecast_interp(
        df_nwp, 
        data["start"], 
        data["stop"],
        data["T_mpc"],
        mpc.Δt,
        k);

    μ_power_L, σ_power_L, μ_wind_L, σ_wind_L = simulate_ode(
        data["T_mpc"], data, itp_mpc_nwp, itp_mpc_dnwp, mpc.Δt, k);

    samples = []
    for i in mpc.𝒩
        push!(samples, lamperti_to_natural_transform_power(μ_power_L[i], data["power_pot_l"][k+i], data["power_pot"][k+i]))
    end

    return(samples)
    
end

function get_simulation_demand(Δt::Int, days::Int, factor::Real)
    hourly_demand_1day = vcat(
        repeat([3.5], 3)..., # demand profile in MW, 3.5 MW from 00:00 to 05:00
        [5.0, 5.25, 6.0, 7.0, 7.5, 6.0]..., # demand profile in MW, 4 MW from 05:00 to 11:00
        [6.0, 5.5, 5.0, 4.5]..., # demand profile in MW, 6 MW from 11:00 to 13:00
        [5.0, 4.5, 6.0, 6.5, 7.0, 7.5, 7.5, 7.0, 6.0]..., # demand profile in MW, 3.5 MW from 13:00 to 22:00
        [5.0 4.5]... # demand profile in MW, 3.5 MW from 22:00 to 00:00
    )

    hourly_demand_1day = vcat(fill.(hourly_demand_1day, Int(60/Δt))...)
    hourly_demand_1day_normalized = hourly_demand_1day/maximum(hourly_demand_1day)
    hourly_demand_1day = hourly_demand_1day_normalized * factor
    return repeat(hourly_demand_1day, days)
end

function get_spot_data(start::DateTime, stop::DateTime, Δt::Int; plotting::Bool = false)
    path="Simulator/data/spot_NO5.csv"
    data = CSV.read(path, DataFrame)
    data[!, "referenceTime"] = DateTime.(data[!, "referenceTime"],DateFormat("yyyy-mm-dd H:M:S"))
    k1 = findfirst(isequal(start), data.referenceTime)
    k2 = findfirst(isequal(stop), data.referenceTime)
    N = k2-k1

    dt = getfield.(Second.(data.referenceTime[k1:k2].-start), :value)
    itp = interpolate((dt,) , data.value[k1:k2], Gridded(Constant(Previous)))

    dt = getfield.(Second.(collect(start:Second(Δt):stop) .- start), :value)
    spot = itp(dt)
    
    if plotting
        p = plot(spot)
        display(p)
    end

    spot_data = Dict(
        "buy" => spot * 3,
        "sell" => spot
    )
    return spot_data
end