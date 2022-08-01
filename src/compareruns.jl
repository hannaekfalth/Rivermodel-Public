using Colors, StatsPlots, CategoricalArrays, Plots.PlotMeasures

# 2000-2019 exist for inflow, 2011-2019 for spot prices, 2004-2019 for production data
function compareruns(years; runfamily="", prod=false, prodplant=false, prodturb=false, choose_max_AandB=true,
                    PDC=false, prodprice=false, prodpriceall=false, taylor=false, recalcprofit=false,
                    plantinv=false, metrics=false, recalculated=false, showplants=false, plantnum=1, runletters=[], plants=[])
    
    params = Dict(y => read_inputdata(y, :standard, silent=true) for y in years)
    
    # nice colors: 
    # :firebrick3, :darkorange1, :gold, :olivedrab, :yellowgreen, 
    # :cyan3, :darkcyan, :darkslateblue, :purple, :ivory4
    
    plotlayouts_list = [   #name, color, markershape, linestyle, fillalpha, runorder         
        "type=NLP, power=bilinear HeadE, e=ncv poly rampseg, head=standardbasic, start=(type=LP, power=E taylor, e=cv segments origo, head=standardbasic)" => ("A_B:L", :darkblue, :circle, :solid, 1.0)
        "type=NLP, power=bilinear HeadE, e=ncv poly rampseg, head=standardbasic, start=(type=LP, power=E constant head, e=cv segments origo, head=constantmax)" => ("A_Cmax", :blue, :circle, :solid, 1.0)
        "type=NLP, power=bilinear HeadE, e=ncv poly rampseg, head=standardbasic, start=(type=LP, power=E constant head, e=cv segments origo, head=constantmean)" => ("A_Cmean", :blue, :circle, :solid, 1.0)
        "type=NLP, power=bilinear HeadE, e=ncv poly rampseg, head=standardbasic, start=(type=MIP, power=bilinear HeadE, e=ncv segments zeroseg, head=standardbasic)" => ("A_MIP", :lightblue, :circle, :solid, 1.0)
        "type=NLP, power=bilinear HeadE, e=cv poly origo, head=standardbasic, start=(type=LP, power=E taylor, e=cv segments origo, head=standardbasic)" => ("B_B:L", :yellowgreen, :circle, :solid, 1.0) 
        "type=NLP, power=bilinear HeadE, e=cv poly origo, head=standardbasic, start=(type=LP, power=E constant head, e=cv segments origo, head=constantmax)" => ("B_Cmax", :green, :circle, :solid, 1.0) 
        "type=NLP, power=bilinear HeadE, e=cv poly origo, head=standardbasic, start=(type=LP, power=E constant head, e=cv segments origo, head=constantmean)" => ("B_Cmean", :green, :circle, :solid, 1.0) 
        "type=LP, power=E taylor, e=cv segments origo, head=standardbasic" => ("B:L", :yellowgreen, :square, :dot, 0.5)
        "type=LP, power=E constant head, e=cv segments origo, head=constantmean" => ("C", :gold, :circle, :solid, 1.0)
        "type=LP, power=E constant head, e=constant eta, head=constantmean" => ("D", :darkorange1, :circle, :solid, 1.0)
        "type=LP, power=aggregated, e=cv segments origo, head=standardbasic" => ("E", :firebrick3, :circle, :solid, 1.0)
        "Historical Data" => ("Historical Data", :grey30, :circle, :solid, 1.0)
    ]
    plotlayouts = Dict(plotlayouts_list)
    runs = readresults(runfamily, years, params, plotlayouts)

    if choose_max_AandB
        A_ = [k for (k,v) in plotlayouts if contains(v[1], "A_")]
        B_ = [k for (k,v) in plotlayouts if contains(v[1], "B_")]
    
        for y in years
            A_profit = Dict()
            B_profit = Dict()
            for run in runs
                if run.name in A_ && run.year==y
                    A_profit[run.name] = run.profit
                elseif run.name in B_
                    B_profit[run.name] = run.profit
                end
            end
            a = findmax(A_profit)[2] 
            b = findmax(B_profit)[2]

            for (i,run) in enumerate(runs)
                if run.name == a && run.year==y
                    runs[i] = (runs[i]..., name="bestArun")
                elseif run.name == b && run.year==y
                    runs[i] = (runs[i]..., name="bestBrun")
                end
            end
        end 

        for n in vcat(A_,B_)
            filter!(r -> !(r.name == n), runs)
            delete!(plotlayouts, n)
        end

        plotlayouts["bestArun"] = ("A", :dodgerblue2, :circle, :solid, 1.0)
        plotlayouts["bestBrun"] = ("B", :yellowgreen, :circle, :solid, 1.0)
    end

    for (key, val) in plotlayouts
        (val[1] == "A" || val[1] == "E" || val[1] == "Historical Data" || endswith(key, "recalculated")) && continue
        plotlayouts["$key, recalculated"] = ("$(val[1])(r)", val[2:4]..., val[5]/2)
    end

    sort!(runs, by = r -> (r.year, get(plotlayouts, r.name, ("Z",))[1]))

    if recalculated || recalcprofit
        recalcargs = (type=:NLP, power="bilinear HeadE", e="ncv poly rampseg", head=:standardbasic)
        recalcstring = choose_max_AandB ? "bestArun" : "type=NLP, power=bilinear HeadE, e=ncv poly rampseg, head=standardbasic"
        newruns = eltype(runs)[]
        for run in runs
            push!(newruns, run)
            if !startswith(run.name, recalcstring) && !startswith(run.name, "type=LP, power=aggregated")
                @unpack name, year, args, TIME, profit, head, power, level, discharge, power_turbines, discharge_turbines, warmstart, solvetime = run
                @unpack PLANT, TURBINE, plant = params[year]
                Water_level = @aa level, TIME, PLANT, [:forebay, :tail]
                Head = @aa head, TIME, PLANT
                results = (Discharge=discharge_turbines, Head, Power_production=power_turbines, Eff_discharge=power_turbines, Water_level)  # only size of Eff_discharge is read by recalculate_variables
                Eff_discharge, power_turbines, profit, level, head = recalculate_variables(params[year], results; recalcargs..., isaggregated=false)
                power = [plant[p].capacity == 0 ? 0.0 : sum(power_turbines[t,p,j] for j in TURBINE[p]) for t in TIME, p in PLANT]
                newrun = (name = "$name, recalculated", year, args, TIME, profit, head, power, level, discharge, power_turbines, discharge_turbines, warmstart, solvetime)
                push!(newruns, newrun)
            end
        end
        runs = newruns
    end

    plotly()
    prod && plot_production(years, runs, params, plotlayouts)
    #prod && plot_production_comparetoA(years, runs, params, plotlayouts)
    prodplant && plot_production_plant(years, runs, params, plotlayouts; runletters, plants)
    prodturb && plot_production_turbine(years, runs, params, plotlayouts; runletters, plants)
    PDC && plot_power_duration(years, runs, params, plotlayouts)
    PDC && showplants && plot_power_duration_plant(years, runs, params, plotlayouts; plants)
    prodprice && !showplants && plot_prod_price(years, runs, params, plotlayouts)
    prodprice && showplants && plot_prod_price_plant(years, runs, params, plotlayouts; plants)
    prodprice && plot_prod_price_article(years, runs, params, plotlayouts)
    prodpriceall && plot_prod_price_all(years, runs, params, plotlayouts)
    # plantinv && plot_plant_investigation(years, runs, params, plotlayouts)
    metrics && metricsvsCT(years, runs, params, plotlayouts)
    metrics && articlemetrics(years, runs, params, plotlayouts) 
    metrics && multiple_run_metrics(years, runs, params, plotlayouts)
    taylor && taylor_eval(years, runs, params, plotlayouts)
    recalcprofit && plotprofit(years, runs, params, plotlayouts)
    
    nothing
end

function listruns(runfamily, years)
    files = readdir("$DATAFOLDER/results")
    matches = match.(Regex("hydrorun (.+?) (\\d+) $(runfamily) ?(type.+?).jld2"), files)
    runs = [(filename=x.match, name=x[3], runtype=x[1], year, args=parse_arguments(x[3])) for x in matches[matches .!== nothing]
                for year = parse(Int, x[2]) if year in years]
    return runs     # Vector of NamedTuples (filename, runtype, year, args)
end

function parse_arguments(argstring)
    startmatch = match(r"(.*?), start=\((.+)\)", argstring)
    if startmatch === nothing
        main, start = argstring, (;)
    else
        main = startmatch[1]
        m = match(r"type=(.+?), power=(.+?), e=(.+?), head=(.+?)$", startmatch[2])
        start = (type=Symbol(m[1]), power=m[2], e=m[3], head=m[4])
    end
    m = match(r"type=(.+?), power=(.+?), e=(.+?), head=(.+?)$", main)
    type, power, e, head = (type=Symbol(m[1]), power=m[2], e=m[3], head=m[4])
    return (; type, power, e, head, start)
end

function get_warmstart_profit(args, year, params)
    type, pow, e, head = (; args..., args.start...)
    recalcargs = (type=args.type, power=args.power, e=args.e, head=args.head)
    argstring = "type=$type, power=$pow, e=$e, head=$head"
    startrun = "$DATAFOLDER/results/hydrorun start $year $argstring.jld2"
    if isfile(startrun)
        head, power_turbines, level, profit, discharge_turbines, solvetime =
                load(startrun, "head", "power_production", "water_level", "profit", "discharge", "solve_t")
        @unpack TIME, PLANT, TURBINE, plant = params[year]
        Water_level = @aa level, TIME, PLANT, [:forebay, :tail]
        Head = @aa head, TIME, PLANT
        results = (Discharge=discharge_turbines, Head, Power_production=power_turbines, Eff_discharge=power_turbines, Water_level)  # only size of Eff_discharge is read by recalculate_variables
        Eff_discharge, power_turbines, profit, level, head = recalculate_variables(params[year], results; recalcargs..., isaggregated=false)
    else
        profit = nothing
    end
    return profit
end

# hardcoded values for now
function MIP_upper_bounds(letter, year)
    !(letter in ["A", "B"] && year in 2015:2019) && return NaN
    row = letter == "A" ? 1 : 2
    col = year-2014
    bounds = [  # 2015:2019
        1176.7   1616.2   1606.8    2353.3    1861.1  # "A"
        1178.6   1618.3   1610.7    2359.3    1873.9  # "B"
    ]
    return bounds[row,col]
end

function readresults(runfamily, years, params, plotlayouts)
    @unpack PLANT, TURBINE, plant = params[years[1]]
    runs = listruns(runfamily, years)
    folder = "$DATAFOLDER/results"
    rundata = NamedTuple[]
    plants = [1:2; 4:16]
    visited = Dict()

    for (filename, name, runtype, year, args) in runs
        !(name in keys(plotlayouts)) && continue
        haskey(visited, (name,year)) && continue
        if args.power == "aggregated"
            power, profit, solvetime = load("$folder/$filename", "power_production", "profit", "solve_t")
            head, level, discharge, power_turbines, discharge_turbines, warmstart = nothing, nothing, nothing, nothing, nothing, nothing
            hours = extrema(key[1] for key in keys(power))  # get TIME from power results
            TIME = hours[1]:hours[2]
        else
            head, power_turbines, level, profit, discharge_turbines, solvetime =
                    load("$folder/$filename", "head", "power_production", "water_level", "profit", "discharge", "solve_t")
            hours = extrema(key[1] for key in keys(power_turbines))  # get TIME from power results
            TIME = hours[1]:hours[2]
            power = AxisArray(0.0, TIME, PLANT)
            discharge = AxisArray(0.0, TIME, PLANT)
            power[:, plants] = [sum(power_turbines[t,p,j] for j in TURBINE[p]) for t in TIME, p in PLANT[plants]]
            discharge[:, plants] = [sum(discharge_turbines[t,p,j] for j in TURBINE[p]) for t in TIME, p in PLANT[plants]]
            warmstart = contains(name, "start") ? get_warmstart_profit(args, year, params) : nothing
        end
        push!(rundata, (; name, year, args, TIME, profit, head, power, level, discharge, power_turbines, discharge_turbines, warmstart, solvetime))
        visited[name, year] = true
    end
    
    return rundata  
end
