using StatsBase

function savevariables(model_year, results, runtype, runfamily, solvetime; type, power, e, head, start=(;))
    println("\nSaving variables to JLD2 archive...")

    folder = "$DATAFOLDER/results"
    runfamily = isempty(runfamily) ? "" : "$(runfamily) "
    basevars = (
        spillage = value.(results.Spillage).data,
        power_production = value.(results.Power_production).data,
        profit = value(results.Profit),
        solve_t = solvetime
    )
    vars = (power == "aggregated") ? basevars : (
        basevars...,
        discharge = value.(results.Discharge).data,
        head = value.(results.Head).data,
        water_level = value.(results.Water_level).data
    )
    argstring = "type=$type, power=$power, e=$e, head=$head"
    startstring = isempty(start) ? "" : ", start=(type=$(start.type), power=$(start.power), e=$(start.e), head=$(start.head))"
    jldsave("$folder/hydrorun $runtype $model_year $runfamily$argstring$startstring.jld2"; vars..., compress=true)
end

function printbasicresults(params, results; type, power, e, head, isaggregated, recalculate=false)
    @unpack rivermodel, model_year, Power_production = results
    println("\nYear: ", model_year)
    println("Profit (Million SEK): ", objective_value(rivermodel))
    try
        println("Best bound: ", objective_bound(rivermodel))
    catch
        println("Best bound: [none available]")
    end
    println("Total power production (TWh): ", sum(value.(Power_production))/1000000)
    if recalculate
        _, powerproduction, profit, _, _ = recalculate_variables(params, results; type, power, e, head, isaggregated)
        println("Recalculated profit (Million SEK): ", profit)
        println("Recalculated power production (TWh): ", sum(powerproduction)/1000000)
    end
end

function printmetrics(params, results)
    @unpack PLANT, plant, hist_prod = params

    power = value.(results.Power_production)
    hourly_total_power = ndims(power) > 1 ? sumdrop(power.data, dims=2) : power.data
    annual_total_power = sum(hourly_total_power)/1000000   #TWh
    annual_historical_production = sum(hist_prod)/1000000 #TWh
    hourly_historical_production = sumdrop(hist_prod, dims=2) #MWh 
    ramping = maximum(abs.(diff(hourly_total_power)))
    ramping_hist = maximum(abs.(diff(hourly_historical_production)))
    
    installed_capacity = [plant[p].capacity for p in PLANT]

    println("\nMetrics: ")
    println("Total power production - Model (TWh): ", annual_total_power)
    println("Total power production - Historical Data (TWh): ", annual_historical_production)
    println()
    println("Maximum Ramping - Model (MWh/h): ", ramping)
    println("Maximum Ramping - Historical Data (MWh/h): ", ramping_hist)
    println()
    println("Maximum Production - Model (MW): ", maximum(hourly_total_power))
    println("Maximum Production - Historical Data (MW): ", maximum(hourly_historical_production))
    println("Minimum Production - Model (MW): ", minimum(hourly_total_power))
    println("Minimum Production - Historical Data (MW): ", minimum(hourly_historical_production))
    println()
    println("Prod. above 85% of total capacity - Model (nr of h): ",
        count(x->(x>= 0.85*sum(installed_capacity)), hourly_total_power)) 
    println("Prod. above 85% of total capacity - Historical Data (nr of h): ",
        count(x->(x>= 0.85*sum(installed_capacity)), hourly_historical_production)) 
    println("Prod. below 15% of total capacity - Model (nr of h): ",
        count(x->(x<= 0.15*sum(installed_capacity)), hourly_total_power)) 
    println("Prod. below 15% of total capacity - Historical Data (nr of h): ",
        count(x->(x<= 0.15*sum(installed_capacity)), hourly_historical_production)) 
    println("No producton - Model (nr of h): ",
        count(x->(x<= 0.0), hourly_total_power)) 
    println("No production - Historical Data (nr of h): ",
        count(x->(x<= 0.0), hourly_historical_production)) 
    println()
    println("RMSD, hourly values: ",
        rmsd(hourly_historical_production, hourly_total_power; normalize=true))
    println("RMSD, weekly moving average: ",
        rmsd(moving_average(hourly_historical_production, 7), moving_average(hourly_total_power, 7); normalize=true))
    println("RMSD, monthly moving average: ",
        rmsd(moving_average(hourly_historical_production, 30), moving_average(hourly_total_power, 30); normalize=true))
end
