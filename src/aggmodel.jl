function buildmodel_agg(params, model_year)
    println("Building model...")
    @unpack PLANT, TIME, plant, aggcap, inflow, spot_price, hist_prod, penalty, dens, grav = params  

    installed_capacity = aggcap
    total_hourly_inflow = sumdrop(inflow, dims=2)
    annual_prod = sum(hist_prod[:,:])
    inflow_MWh = total_hourly_inflow./sum(total_hourly_inflow)*annual_prod
    maximum_reservoir_capacity = 0 #TWh 
    for p in PLANT
        maximum_reservoir_capacity =
            maximum_reservoir_capacity + plant[p].reservoir*60*60*dens*grav*plant[p].reportedhead/1e6*0.9
    end
    
    rivermodel = Model(Gurobi.Optimizer)
    set_optimizer_attributes(rivermodel, "Threads" => 12, "Method" => 2, "ScaleFlag" => -1,
        "NonConvex" => 2, "MIPGap" => 1e-9)#, "OutputFlag" => 0)
    # m = Model(solver=GurobiSolver(Method=2, Threads=threads, BarConvTol=1.5e-8, Crossover=0)) 
    # BarConvTol=1e-8, NumericFocus=2 (0-3), ScaleFlag=3 (-1-3), Crossover=0

    @variables rivermodel begin
    	Profit					   									# Million SEK
        Reservoir_content[t in TIME]      		>=0                 # MWh
        Power_production[t in TIME]       		>=0					# MWh
        Spillage[t in TIME]	       		        >=0					# MWh
    
    end #variables

    for t in TIME
        set_upper_bound(Reservoir_content[t], maximum_reservoir_capacity)
        set_upper_bound(Power_production[t], installed_capacity)
    end

    if penalty > 0
        @variables rivermodel begin
            DischargeDelta[t in TIME]                   >=0    # MWh/h
        end #variables

        @constraints rivermodel begin
        # Only one of these two constraints will be active. 
        Calculate_DischargeDelta1[t in TIME],
            DischargeDelta[t] >= Power_production[t] - Power_production[shift(TIME, t-1)]

        Calculate_DischargeDelta2[t in TIME, p in POWERPLANT],
            DischargeDelta[t] >= Power_production[shift(TIME, t-1)] - Power_production[t]
        end #constraints
    end

    initial_content = 0.70 * maximum_reservoir_capacity # innehåll i reservoaren timme 0

    @constraints rivermodel begin
        Reservoir_final,
            Reservoir_content[length(TIME)] == initial_content
        
        Water_Balance[t in TIME],   
            Reservoir_content[t] == (t>1 ? Reservoir_content[t-1] : initial_content) + 
                + inflow_MWh[t] - Spillage[t] - Power_production[t]	

        Calculate_Profit,
            Profit == sum(Power_production[t]*spot_price[t] +
                - ((penalty > 0) ? penalty * DischargeDelta[t] : 0) for t in TIME)/1000000
   
    end #constraints

    @objective rivermodel Max begin
        		Profit
    end #objective

    return (; rivermodel, model_year, Power_production, Reservoir_content, Spillage, Profit)
end
