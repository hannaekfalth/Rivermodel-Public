function buildmodel(params, model_year; type, power, e, head)
    println("\nBuilding model...")
    println("year = ", model_year)
    println("type = ", type)
    println("power = ", power)
    println("e = ", e)
    println("head = ", head)
    @unpack PLANT, TURBINE, SEGMENT, TIME, POINT, LAGS, plant, downstream, inflow, etacoeff, spot_price, penalty,
        reservoir_start, reservoir_end, reservoir_area, grav, dens, WtoMW,
        maxturbdischarge, maxturbeta, maxplantdischarge, meanturbdischarge, meanturbeta, maxturbeta, k_firstseg, end_rampseg, end_zeroseg, end_zeroseg_poly, end_origoseg,
        k_segcoeff, m_segcoeff, k_segcoeff_origo, m_segcoeff_origo, bigM,
        minlevel, maxlevel, tailrace_meanlevel, maxhead, meanhead, minhead,
        tailrace_per_dischargelags, tailrace_per_forebay, tailrace_per_downstreamforebay,
        tailrace_constant, hist_prod, hist_discharge, hist_head, hist_spill, hist_level, minspill, minflow, vars, plants = params

    hoursperperiod = 1
    
    rivermodel = Model()
    
    @variables rivermodel begin
        Profit                                                      # unbounded, it's the objective                                     # million SEK
        Power_production[t in TIME, p in PLANT, j in TURBINE[p]]    # bounds set below                                                  # MWh/period
        Reservoir_content[t in TIME, p in PLANT]                    >= 0,                       (upper_bound = plant[p].reservoir)      # HE (m3/s * 1h)
        Water_level[t in TIME, p in PLANT, i in POINT]              >= minlevel[t,p,i],         (upper_bound = maxlevel[t,p,i])         # m  
        Discharge[t in TIME, p in PLANT, j in TURBINE[p]]           >= 0,                       (upper_bound = maxturbdischarge[p,j])   # m3/s
        Spillage[t in TIME, p in PLANT]                             >= minspill[t,p],           (upper_bound = 600)                     # m3/s
        Head[t in TIME, p in PLANT]                                 >= minhead[p],              (upper_bound = maxhead[p])              # m 
        Eff_discharge[t in TIME, p in PLANT, j in TURBINE[p]]                  >= 0,                       (upper_bound = maxturbdischarge[p,j])
    end #variables

    if type == :MIP
        @variable(rivermodel, Z[t in TIME, p in PLANT, j in TURBINE[p]], Bin)
    end

    if !contains(power, "taylor")
        for t in TIME, p in PLANT, j in TURBINE[p]
            set_lower_bound(Power_production[t,p,j], 0)
        end
        # should add upper bound for all runs
    end

    @constraints rivermodel begin
        Reservoir_final[p in PLANT],
        	Reservoir_content[TIME[end], p] == reservoir_end[p]	
        
        Water_Balance[t in TIME, p in PLANT],
            Reservoir_content[t,p] == (t > TIME[1] ? Reservoir_content[t-1, p] : reservoir_start[p]) +
                hoursperperiod * (
                    + inflow[t,p] - Spillage[t,p] - sum(Discharge[t,p,j] for j in TURBINE[p])  + 
                    + sum(up.flowshare * sum(Discharge[shift(TIME, t - up.dischargedelay), up.name, j] for j in TURBINE[up.name])
                                for up in plant[p].upstream) + 
                    + sum(up.flowshare * Spillage[shift(TIME, t - up.spilldelay), up.name]
                                for up in plant[p].upstream)
                )
        
        Forebay_level[t in TIME, p in PLANT],
            Water_level[t, p, :forebay] == plant[p].reservoirlow + Reservoir_content[t, p] / reservoir_area[p]          #This gives us an upper and lower bound on forebay water level
        
        Calculate_head[t in TIME, p in PLANT],
            Head[t,p] == Water_level[t, p, :forebay] - Water_level[t, p, :tail]

        Minimum_flow[t in TIME, p in PLANT],	  
        	sum(Discharge[t,p,j] for j in TURBINE[p]) + Spillage[t,p] >= minflow[t,p]

        Calculate_Profit,
            Profit == sum(sum(Power_production[t,p,j] for p in PLANT for j in TURBINE[p]) * spot_price[t] for t in TIME) / 1e6
    end #constraints
    
    if head == :regression || startswith(string(head), "standard")
        @constraints rivermodel begin
            Tailrace_level[t in TIME, p in PLANT],
                Water_level[t, p, :tail] == tailrace_constant[p] +
                    sum(tailrace_per_dischargelags[p,n+1] * sum(Discharge[shift(TIME, t-n), p, j] for j in TURBINE[p]) for n in LAGS) +
                    + tailrace_per_forebay[p] * Water_level[t, p, :forebay] +
                        ((p == PLANT[end]) ? 0.0 :
                            tailrace_per_downstreamforebay[p] * Water_level[t, downstream[p], :forebay])
        end
    elseif head == :singlereservoir
        @constraints rivermodel begin
            Tailrace_level[t in TIME, p in PLANT],
                Water_level[t, p, :tail] == tailrace_meanlevel[p]  
        end
    elseif head == :reservoirdiff
        @constraints rivermodel begin
            Tailrace_level[t in TIME, p in PLANT],
                # this assumes there is only one downstream plant (i.e. no forks)
                Water_level[t, p, :tail] ==
                    ((p == PLANT[end]) ? 0 : Water_level[t, downstream[p], :forebay])
        end
    end

    if e == "cv segments origo"
        @constraints rivermodel begin
            eta_segments[t in TIME, p in PLANT, j in TURBINE[p], s in SEGMENT], 
                Eff_discharge[t,p,j] <= k_segcoeff_origo[p,j,s] * Discharge[t,p,j] + m_segcoeff_origo[p,j,s]
        end
    elseif e == "cv segments noseg"
        @constraints rivermodel begin
            eta_segments[t in TIME, p in PLANT, j in TURBINE[p], s in SEGMENT], 
                Eff_discharge[t,p,j] <= k_segcoeff[p,j,s] * Discharge[t,p,j] + m_segcoeff[p,j,s]
        end
    elseif e == "constant eta"
        @constraints rivermodel begin
            eta_constant[t in TIME, p in PLANT, j in TURBINE[p]],
                Eff_discharge[t,p,j] <= Discharge[t,p,j] * maxturbeta[p,j]
        end
    elseif type == :NLP && e == "ncv segments rampseg"
        @NLconstraint(rivermodel,
            eta_segments[t in TIME, p in PLANT, j in TURBINE[p], s in SEGMENT], 
                Eff_discharge[t,p,j] <= (Discharge[t,p,j] <= end_rampseg[p,j]) * k_firstseg[p,j] * Discharge[t,p,j] +
                    (Discharge[t,p,j] > end_rampseg[p,j]) * (k_segcoeff[p,j,s] * Discharge[t,p,j] + m_segcoeff[p,j,s])
        )
    elseif type == :MIP && e == "ncv segments rampseg"
        @constraints rivermodel begin
            eta_segments[t in TIME, p in PLANT, j in TURBINE[p], s in SEGMENT], 
                Eff_discharge[t,p,j] <= k_segcoeff[p,j,s] * Discharge[t,p,j] + m_segcoeff[p,j,s] + bigM[p,j] * (1 - Z[t,p,j])

            eta_segments_ramp[t in TIME, p in PLANT, j in TURBINE[p]], 
                Eff_discharge[t,p,j] <= k_firstseg[p,j] * Discharge[t,p,j] + bigM[p,j] * Z[t,p,j]

            eta_segments_zero[t in TIME, p in PLANT, j in TURBINE[p]], 
                Discharge[t,p,j] <= end_rampseg[p,j] + (maxturbdischarge[p,j] - end_rampseg[p,j]) * Z[t,p,j]
        end
    elseif type == :MIP && e == "ncv segments zeroseg"
        @constraints rivermodel begin
            eta_segments[t in TIME, p in PLANT, j in TURBINE[p], s in SEGMENT], 
                Eff_discharge[t,p,j] <= k_segcoeff[p,j,s] * Discharge[t,p,j] + m_segcoeff[p,j,s] + bigM[p,j] * (1 - Z[t,p,j])

            eta_segments_ramp[t in TIME, p in PLANT, j in TURBINE[p]], 
                Eff_discharge[t,p,j] <= bigM[p,j] * Z[t,p,j]

            eta_segments_zero[t in TIME, p in PLANT, j in TURBINE[p]], 
               Discharge[t,p,j] <= maxturbdischarge[p,j] * Z[t,p,j]

            eta_segments_power[t in TIME, p in PLANT, j in TURBINE[p]],
               Power_production[t,p,j] <= maxhead[p] * bigM[p,j] * Z[t,p,j] * grav * dens * WtoMW
        end
    elseif type == :NLP && e == "cv poly origo"
        @NLconstraint(rivermodel,
            eta_poly[t in TIME, p in PLANT, j in TURBINE[p]],
                Eff_discharge[t,p,j] <= (Discharge[t,p,j] <= end_origoseg[p,j]) * k_segcoeff_origo[p,j,1] * Discharge[t,p,j] +
                    (Discharge[t,p,j] > end_origoseg[p,j]) * sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
        )
    elseif type == :MIP && e == "cv poly noseg"
        @constraint(rivermodel,
            eta_poly[t in TIME, p in PLANT, j in TURBINE[p]], 
                Eff_discharge[t,p,j] <= sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
        )
    elseif type == :NLP && e == "cv poly noseg"
        @NLconstraint(rivermodel,
            eta_poly[t in TIME, p in PLANT, j in TURBINE[p]], 
                Eff_discharge[t,p,j] <= sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
        )
    elseif type == :NLP && e == "ncv poly rampseg"
        @NLconstraint(rivermodel,
            eta_poly[t in TIME, p in PLANT, j in TURBINE[p]],
                Eff_discharge[t,p,j] <= (Discharge[t,p,j] <= end_rampseg[p,j]) * k_firstseg[p,j] * Discharge[t,p,j] +
                    (Discharge[t,p,j] > end_rampseg[p,j]) * sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
        )
    elseif type == :NLP && e == "ncv poly zeroseg"
        @NLconstraint(rivermodel,
            eta_poly[t in TIME, p in PLANT, j in TURBINE[p]],
                Eff_discharge[t,p,j] <= (Discharge[t,p,j] > end_zeroseg_poly[p,j]) * sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
        )
    else 
        error("No alternative with: type = $type, e = $e.")
    end
    
    if head == :constantmean || head == :constantmax
        @constraints rivermodel begin
            Calc_Power_production[t in TIME, p in PLANT, j in TURBINE[p]],
                Power_production[t,p,j] == max(0, head == :constantmean ? meanhead[p] : maxhead[p]) * Eff_discharge[t,p,j] * grav * dens * WtoMW
        end
    elseif contains(power, "taylor")
        @constraints rivermodel begin
            Calc_Power_production[t in TIME, p in PLANT, j in TURBINE[p]],
                Power_production[t,p,j] == (plant[p].capacity == 0 ? 0.0 :
                    (meanhead[p] * Eff_discharge[t,p,j] + Head[t,p] * meanturbdischarge[p,j] * meanturbeta[p,j] +
                        - meanhead[p] * meanturbeta[p,j] * meanturbdischarge[p,j]) * grav * dens * WtoMW)
        end
    elseif type == :NLP && power == "bilinear HeadE"
        @NLconstraint(rivermodel,
            [t in TIME, p in PLANT, j in TURBINE[p]], # Workaround for Gurobi bug, not naming constraint saves a lot of model generation time
                Power_production[t,p,j] <= Head[t,p] * Eff_discharge[t,p,j] * grav * dens * WtoMW
        )
    elseif type == :MIP && power == "bilinear HeadE"
        @constraints rivermodel begin
            [t in TIME, p in PLANT, j in TURBINE[p]],   # Workaround for Gurobi bug, not naming constraint saves a lot of model generation time
                Power_production[t,p,j] <= Head[t,p] * Eff_discharge[t,p,j] * grav * dens * WtoMW
        end
    else 
        error("No alternative with: type = $type, power = $power, head = $head.")
    end

    @objective(rivermodel, Max, Profit)

    return (; rivermodel, model_year, Power_production, Head, Discharge, Eff_discharge, Spillage,
                Reservoir_content, Water_level, Profit)    #, Calc_Power_production)
end



# Set start values using equations in run 2, no matter how variables in run 1 were calculated.
function set_start_values!(params, results, results2; type, power, e, head)
    @unpack TIME, PLANT, TURBINE = params
    Eff_discharge, Power_production, Profit, Water_level, Head = recalculate_variables(params, results; type, power, e, head, isaggregated=false)
    for t in TIME, p in PLANT
        set_start_value(results2.Water_level[t,p,:tail], Water_level[t,p,:tail])
        set_start_value(results2.Head[t,p], Head[t,p])
        for j in TURBINE[p]
            set_start_value(results2.Eff_discharge[t,p,j], Eff_discharge[t,p,j])
            set_start_value(results2.Power_production[t,p,j], Power_production[t,p,j])
        end
    end
    set_start_value(results2.Profit, Profit)
end

function recalculate_variables(params, results; type, power, e, head, isaggregated)
    @unpack TIME, PLANT, SEGMENT, plant, TURBINE, spot_price, grav, dens, WtoMW, meanhead, maxhead, etacoeff, meanturbeta, maxturbeta, meanturbdischarge,
            end_rampseg, end_zeroseg, end_zeroseg_poly, k_firstseg, k_segcoeff, m_segcoeff, k_segcoeff_origo, m_segcoeff_origo, end_origoseg,
            tailrace_constant, tailrace_per_dischargelags, tailrace_per_forebay, tailrace_per_downstreamforebay, downstream, LAGS = params
    
    if isaggregated
        Power_production, Profit = value.(results.Power_production), value.(results.Profit)
        return nothing, Power_production, Profit, nothing, nothing
    elseif typeof(results.Power_production) <: JuMP.Containers.SparseAxisArray
        Discharge, Head, Eff_discharge, Power_production, Water_level =
            value.(results.Discharge), value.(results.Head), value.(results.Eff_discharge), value.(results.Power_production), value.(results.Water_level)
    else
        @unpack Discharge, Head, Eff_discharge, Power_production, Water_level = results
    end

    for t in TIME, p in PLANT
        if head == :regression || startswith(string(head), "standard")
            Water_level[t, p, :tail] = tailrace_constant[p] +
                sum(tailrace_per_dischargelags[p,n+1] * (plant[p].capacity == 0 ? 0.0 : sum(Discharge[shift(TIME, t-n), p, j] for j in TURBINE[p])) for n in LAGS) +
                    + tailrace_per_forebay[p] * Water_level[t, p, :forebay] +
                            ((p == PLANT[end]) ? 0.0 :
                                tailrace_per_downstreamforebay[p] * Water_level[t, downstream[p], :forebay])
        elseif head == :singlereservoir
            Water_level[t, p, :tail] = tailrace_meanlevel[p]
        elseif head == :reservoirdiff
            Water_level[t, p, :tail] = ((p == PLANT[end]) ? 0 : Water_level[t, downstream[p], :forebay])
        end

        Head[t,p] = Water_level[t, p, :forebay] - Water_level[t, p, :tail]

        for j in TURBINE[p]
            if e == "cv segments origo"
                Eff_discharge[t,p,j] = minimum(k_segcoeff_origo[p,j,s] * Discharge[t,p,j] + m_segcoeff_origo[p,j,s] for s in SEGMENT)
            elseif e == "cv segments noseg"
                Eff_discharge[t,p,j] = minimum(k_segcoeff[p,j,s] * Discharge[t,p,j] + m_segcoeff[p,j,s] for s in SEGMENT)
            elseif e == "constant eta"
                Eff_discharge[t,p,j] = Discharge[t,p,j] * maxturbeta[p,j]
            elseif e == "ncv segments rampseg"
                Eff_discharge[t,p,j] = (Discharge[t,p,j] <= end_rampseg[p,j]) ? k_firstseg[p,j] * Discharge[t,p,j] :
                        minimum(k_segcoeff[p,j,s] * Discharge[t,p,j] + m_segcoeff[p,j,s] for s in SEGMENT)
            elseif e == "ncv segments zeroseg"
                Eff_discharge[t,p,j] = (Discharge[t,p,j] <= end_zeroseg[p,j]) ? 0.0 :
                        minimum(k_segcoeff[p,j,s] * Discharge[t,p,j] + m_segcoeff[p,j,s] for s in SEGMENT)
            elseif (type == :NLP || type == :MIP) && e == "cv poly noseg"
                Eff_discharge[t,p,j] = sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
            elseif type == :NLP && e == "ncv poly rampseg"
                Eff_discharge[t,p,j] = (Discharge[t,p,j] <= end_rampseg[p,j]) ? k_firstseg[p,j] * Discharge[t,p,j] :
                        sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
            elseif type == :NLP && e == "ncv poly zeroseg"
                Eff_discharge[t,p,j] = (Discharge[t,p,j] <= end_zeroseg_poly[p,j]) ? 0.0 :
                        sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
            elseif type == :NLP && e == "cv poly origo"
                Eff_discharge[t,p,j] = (Discharge[t,p,j] <= end_origoseg[p,j]) * k_segcoeff_origo[p,j,1] * Discharge[t,p,j] +
                    (Discharge[t,p,j] > end_origoseg[p,j]) * sum(etacoeff[p,j][i] * Discharge[t,p,j]^(i-1) for i = 1:3)
            end

            if head == :constantmean || head == :constantmax
                Power_production[t,p,j] = max(0, head == :constantmean ? meanhead[p] : maxhead[p]) * Eff_discharge[t,p,j] * grav * dens * WtoMW
            elseif contains(power, "taylor")
                Power_production[t,p,j] = (plant[p].capacity == 0) ? 0.0 :
                    (meanhead[p] * Eff_discharge[t,p,j] + Head[t,p] * meanturbdischarge[p,j] * meanturbeta[p,j] +
                        - meanhead[p] * meanturbeta[p,j] * meanturbdischarge[p,j]) * grav * dens * WtoMW
            elseif power == "bilinear HeadE"
                Power_production[t,p,j] = Head[t,p] * Eff_discharge[t,p,j] * grav * dens * WtoMW
            end
        end
    end

    Profit = sum(sum(Power_production[t,p,j] for p in PLANT for j in TURBINE[p]) * spot_price[t] for t in TIME) / 1e6

    return Eff_discharge, Power_production, Profit, Water_level, Head
end
