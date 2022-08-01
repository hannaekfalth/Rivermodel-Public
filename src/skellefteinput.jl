using Dates, JLD2, Polynomials

function read_inputdata(model_year, headtype; weeks=1:52, env_con=true, scaleeta=false, etaquantiles=nothing,
                silent=false, regressionyears=2015:2020, regressionvars=Bool[1,1,1,1,1,1,1])
    !silent && println("Reading input data...")

    plants, TURBINE, etapoints, set_environmental_constraints! = riverdata(:Skellefte)
    nplants = length(plants)
    plant = Dict(p.name => p for p in plants)

    TIME = 1:daysinyear(model_year)*24  # overwritten at bottom!!!
    PLANT = [p.name for p in plants]
    POINT = [:forebay, :tail]
    SEGMENT=1:10
    LAGS = 0:3
    downstream = Dict(up.name => p.name for p in plants for up in p.upstream)
    
    vars = load("$DATAFOLDER/skelleftedata_raw.jld2")
    @unpack hours, production, forebay_level, tail_level, discharge, spillage, reservoir_weekly, inflow_daily, elspot_price = vars
    year = Dates.year.(hours)
    inflow_year = Dates.year.(Date("2000-01-01"):Day(1):Date("2019-12-31"))
    inflow = @aa repeat(inflow_daily[inflow_year .== model_year, :], inner = (24,1)), TIME, PLANT
    
    grav = 9.81 # m/s2
    dens = 998 # kg/m3
    WtoMW = 1/1e6
   
    penalty = 0.0  # Makes it possible to put a cost on changing discharge from one hour to the next

    yr = (year .== model_year)

    kvist = findfirst(PLANT .== :Kvistforsen)
    up = @view forebay_level[yr, kvist]
    quick_clean!(up, mean(up[up.>0]))
    low = @view tail_level[yr, kvist]
    quick_clean!(low, mean(low[low.>0]))

    hist_head = @aa forebay_level[yr,:] - tail_level[yr,:], TIME, PLANT
    hist_prod = @aa production[yr,:], TIME, PLANT
    hist_spill = @aa spillage[yr,:], TIME, PLANT
    hist_discharge = @aa discharge[yr,:], TIME, PLANT
    hist_level = @aa 0.0, TIME, PLANT, POINT
    hist_level[:, :, :forebay] = forebay_level[yr,:]
    hist_level[:, :, :tail] = tail_level[yr,:]
    spot_price = elspot_price[yr]

    minlevel = @aa zeros(length(TIME), nplants, length(POINT)), TIME, PLANT, POINT  # m.a.s.l. upper bound on Water_level
    maxlevel = @aa zeros(length(TIME), nplants, length(POINT)), TIME, PLANT, POINT  # m.a.s.l. lower bound on Water_level
    for p in plants
        maxlevel[:, p.name, :forebay] .= p.reservoirhigh
        minlevel[:, p.name, :forebay] .= p.reservoirlow
    end

    if headtype == :reservoirdiff
        for p in PLANT
            maxlevel[:, p, :tail] .= ((p == PLANT[end]) ? 0 : maxlevel[:, downstream[p], :forebay])
            minlevel[:, p, :tail] .= ((p == PLANT[end]) ? 0 : minlevel[:, downstream[p], :forebay])
        end
    else
        for (i,p) in enumerate(PLANT)
            plants[i].capacity == 0 && continue
            maxlevel[:, p, :tail] .= positivequantile(tail_level[:,i], 0.99)    # max m.a.s.l. tailrace
            minlevel[:, p, :tail] .= positivequantile(tail_level[:,i], 0.01)    # min m.a.s.l. tailrace
        end
        maxlevel[:, :Hornavan, :tail] = maxlevel[:, downstream[:Hornavan], :forebay]
        maxlevel[:, :Bergsby, :tail] .= 0
        minlevel[:, :Hornavan, :tail] = minlevel[:, downstream[:Hornavan], :forebay]
        minlevel[:, :Bergsby, :tail] .= 0
    end

#=     if !silent
        println("\nWater level bounds [max forebay, min forebay, max tail, min tail]")
        display([maxlevel[1,:,:forebay] minlevel[1,:,:forebay] maxlevel[1,:,:tail] minlevel[1,:,:tail]])
    end =#
    
    tailrace_meanlevel = @aa [mean(tail_level[tail_level[:,i].>0, i]) for i = 1:nplants], PLANT # m.a.s
    tailrace_meanlevel[:Hornavan] = mean((maxlevel[t, :Hornavan, :tail] + minlevel[t, :Hornavan, :tail])/2 for t in TIME)
    tailrace_meanlevel[:Bergsby] = mean((maxlevel[t, :Bergsby, :tail] + minlevel[t, :Bergsby, :tail])/2 for t in TIME)

    head = forebay_level .- tail_level                                              # m.a.s
    head[(forebay_level .<= 0) .| (tail_level .<= 0)] .= -Inf                       # remove unphysical head values
    maxhead_data = @aa [positivequantile(head[:, i], 0.99) for i = 1:nplants], PLANT     # m.a.s
    minhead_data = @aa [positivequantile(head[:, i], 0.01) for i = 1:nplants], PLANT     # m.a.s
    meanhead = @aa [positivequantile(head[yr, i], 0.5) for i = 1:nplants], PLANT     # m.a.s
    
    tailrace_constant = @aa 0.0, PLANT
    tailrace_per_dischargelags = @aa 0.0, PLANT, LAGS
    tailrace_per_forebay = @aa 0.0, PLANT
    tailrace_per_downstreamforebay = @aa 0.0, PLANT
    meanabserror = zeros(nplants)

    okyears = [y in regressionyears for y in year]
    for i in 1:nplants
        idown = clamp(i+1, 3, nplants)
        regcoeff, meanabserror[i] =
            tailrace_regression(discharge[:,i], forebay_level[:,i], forebay_level[:,idown],
                            tail_level[:,i], okyears, regressionvars, LAGS, i == nplants)
        tailrace_constant[i], tailrace_per_forebay[i], tailrace_per_downstreamforebay[i] =
                regcoeff[[1,6,7]]
        tailrace_per_dischargelags[i, :] = regcoeff[2:5]
    end

    realplants = [plant[p].capacity > 0 for p in PLANT]
#=     if !silent && contains(string(headtype), "standard")
        println("\nTailrace regression errors: mean ", round(mean(meanabserror[realplants]), digits=3))
        display(round.(meanabserror, digits=3)')
    end =#
    
    if startswith(string(headtype), "standard")
        #= if !silent
            println("\nMinimum tail water level, using historical data...")
            # display([minlevel[1,:,:tail] tailrace_constant])
        end =#
        minlevel[:,:Hornavan,:tail] = minlevel[:,downstream[:Hornavan],:forebay]
        minlevel[:,.!realplants,:tail] .= 0 
        #= if !silent
             println("\nNew water level bounds [max forebay, min forebay, max tail, min tail]")
             display([maxlevel[1,:,:forebay] minlevel[1,:,:forebay] maxlevel[1,:,:tail] minlevel[1,:,:tail]])
        end =#
    end

    minflow = @aa zeros(length(TIME), nplants), TIME, PLANT     # m3/s (spill or discharge)
    minspill = @aa zeros(length(TIME), nplants), TIME, PLANT    # m3/s

    if env_con
        !silent && println("\nApplying environmental constraints...")
        set_environmental_constraints!(model_year, minspill, minlevel, minflow)
    end

    for t in TIME[2:end], p in PLANT[1:end-1]
        if hist_level[t,p,:forebay] > maxlevel[t,p,:forebay] 
            hist_level[t,p,:forebay] = hist_level[t-1,p,:forebay]
        elseif hist_level[t,p,:forebay] < minlevel[t,p,:forebay] 
            hist_level[t,p,:forebay] = hist_level[t-1,p,:forebay]
        end
        if hist_level[t,p,:tail] >= hist_level[t,p,:forebay]
            hist_level[t,p,:tail] = hist_level[t-1,p,:tail]
        elseif hist_level[t,p,:tail] < hist_level[t,downstream[p],:forebay]
            hist_level[t,p,:tail] = hist_level[t-1,p,:tail]
        end
    end

    maxhead = @aa [maxlevel[1, p, :forebay] - minlevel[1, p, :tail] for p in PLANT], PLANT  #m.a.s
    minhead = @aa [minlevel[1, p, :forebay] - maxlevel[1, p, :tail] for p in PLANT], PLANT  #m.a.s

    etacoeff, maxturbdischarge, maxturbeta, meanturbdischarge, meanturbeta, maxplantdischarge, 
        k_firstseg, end_rampseg, end_zeroseg, end_zeroseg_poly, end_origoseg, k_segcoeff, m_segcoeff, k_segcoeff_origo, m_segcoeff_origo, bigM =
            estimate_eta_discharge(vars, grav, dens, WtoMW, PLANT, TURBINE, SEGMENT, etapoints, plants, silent, maxhead, minhead, meanhead)
    
    firsthour = 1 + 7*24*(weeks[1] - 1)      # model_year-01-01T00:00:00
    lasthour = weeks[end] < 52 ? 7*24*weeks[end] : daysinyear(model_year)*24

    reservoir_content_hist = @aa [(hist_level[t, p.name, :forebay] - p.reservoirlow)/(p.reservoirhigh - p.reservoirlow)* p.reservoir for t in TIME, p in plants], TIME, PLANT    # HE   
    
        #plotly()
        #for p in PLANT
        #plot(hist_level[:,p,:forebay], title = "$p - forebay level") |> display 
        #plot(reservoir_content_hist[:,p], title = "$p - reservoir content hist") |> display    
        #end 
         
    reservoir_start = @aa [reservoir_content_hist[firsthour,p] for p in PLANT], PLANT   # HE
    reservoir_end = @aa [reservoir_content_hist[lasthour,p] for p in PLANT], PLANT     # HE
    reservoir_area = @aa [p.reservoir / (p.reservoirhigh - p.reservoirlow) for p in plants], PLANT  # HE/m
    
    TIME = firsthour:lasthour

    cap = Dict()
    for p in PLANT
        max_ed = Dict()
        for j in TURBINE[p]
        max_ed[j] = maximum([etapoints[p,j][i].d * etapoints[p,j][i].e for i in 1:length(etapoints[p,j])])
        end
        cap[p] = plant[p].capacity == 0 ? 0.0 : sum(maxhead[p]*max_ed[j]*grav*dens*WtoMW for j in TURBINE[p])
        #println("Maxcap $p = ", cap[p])
    end 
    
    aggcap = sum(cap[p] for p in PLANT)
    #println("Theoretic max capacity in river = ", aggcap)


    return (; PLANT, TURBINE, SEGMENT, TIME, POINT, LAGS, plant, downstream, aggcap, inflow, etacoeff, spot_price, penalty,
        reservoir_start, reservoir_end, reservoir_area, grav, dens, WtoMW,
        maxturbdischarge, maxturbeta, maxplantdischarge, meanturbdischarge, meanturbeta, k_firstseg, end_rampseg, end_zeroseg, end_zeroseg_poly, end_origoseg,
        k_segcoeff, m_segcoeff, k_segcoeff_origo, m_segcoeff_origo, bigM,
        minlevel, maxlevel, tailrace_meanlevel, maxhead, meanhead, minhead,
        tailrace_per_dischargelags, tailrace_per_forebay, tailrace_per_downstreamforebay,
        tailrace_constant, hist_prod, hist_discharge, hist_head, hist_spill, hist_level, minspill, minflow, vars, plants)
end

function tailrace_regression(discharge, forebay, forebay_down, tail, ok, regressionvars, LAGS, lastplant)
    vars = [ones(length(discharge)) lags.(Ref(discharge), LAGS)... forebay forebay_down]
    f = ok .& (tail .> 0) .& all(vars[:, 1:6][:, regressionvars[1:6]] .> 0, dims=2) |> vec
    if lastplant
        vars[:, 7] .= 0.0
    else
        f = f .& (forebay_down .> 0) .& (tail .>= forebay_down)
    end
    vars[:, .!regressionvars] .= 0.0
    regcoeff = vars[f, :] \ tail[f]
    pred_tail = vars[f, :] * regcoeff
    meanabserror = mean(abs.(tail[f] .- pred_tail))
    return regcoeff, meanabserror
end

function estimate_eta_discharge(vars, grav, dens, WtoMW, PLANT, TURBINE, SEGMENT, etapoints, plants, silent, maxhead, minhead, meanhead)
    @unpack hours, production, forebay_level, tail_level, discharge,  = vars

    year = Dates.year.(hours)
    yr = (year .>= 2015) .& (year .<= 2020)

    kvist = findfirst(PLANT .== :Kvistforsen)
    up = @view forebay_level[yr, kvist]
    quick_clean!(up, mean(up[up.>0]))
    low = @view tail_level[yr, kvist]
    quick_clean!(low, mean(low[low.>0]))

    maxplantdischarge = zeros(length(PLANT))
    maxplanteta = zeros(length(PLANT))
    meanplantdischarge = zeros(length(PLANT))
    meanplanteta = zeros(length(PLANT))
    etacoeff = Dict()
    maxturbdischarge = Dict()
    maxturbeta = Dict()
    meanturbdischarge = Dict()
    meanturbeta = Dict()

    inflectionpoint = Dict() 
    origolinepoint = Dict() 
    segmentpoints = Dict()
    origosegmentpoints = Dict()
    k_segcoeff = Dict()
    k_segcoeff_origo = Dict()
    m_segcoeff = Dict()
    m_segcoeff_origo = Dict()
    yiszero = Dict()
    bigM = Dict()
    k_firstseg = Dict()
    end_rampseg = Dict()
    end_zeroseg = Dict()
    end_zeroseg_poly = Dict()
    end_origoseg = Dict()
    Terror = Dict()
    Cerror = Dict()
 
    for (i, p) in enumerate(PLANT), j in TURBINE[p]
        plants[i].capacity == 0 && continue
        
        head = forebay_level[yr, i] - tail_level[yr, i]
        dc0 = discharge[yr, i]
        calc_elec = grav * dens * dc0 .* head * WtoMW
        eta0 = production[yr, i] ./ calc_elec
        f = (eta0 .> 0) .& (eta0 .< 1) .& (forebay_level[yr, i] .> 0) .& (tail_level[yr, i] .> 0) .&
                    (dc0 .> 0) .& (production[yr, i] .> 0)
        dc, eta = dc0[f], eta0[f] # discharge and eta per hour for plant i
        maxplantdischarge[i] = quantile(dc, 0.999)
        maxplanteta[i] = quantile(eta, 0.999)
        meanplantdischarge[i] = quantile(dc, 0.5)
        meanplanteta[i] = quantile(eta, 0.5)
        
        #meanplantdischarge = @aa [positivequantile(discharge[discharge[:,i].>0, i], 0.5) for i = 1:nplants], PLANT

        p_discharge = [etapoints[p,j][i].d for i in 1:length(etapoints[p,j])]
        p_eta = [etapoints[p,j][i].e for i in 1:length(etapoints[p,j])]
        p_pe = p_discharge .* p_eta
        Eff_discharge = Polynomials.fit(p_discharge, p_pe, 2)  # fit 2nd degree polynomial to Eff_discharge points (not eta!)
        etacoeff[p,j] = coeffs(Eff_discharge)
        dPE = derivative(Eff_discharge)
        root = minimum(roots(Eff_discharge))

        maxturbeta[p,j] = maximum(p_eta)
        maxturbdischarge[p,j] = p_discharge[end]*1.05 
        meanturbdischarge[p,j] = meanplantdischarge[i]/maxplantdischarge[i]*maxturbdischarge[p,j] 
        meanturbeta[p,j] = meanplanteta[i]/maxplanteta[i]*maxturbeta[p,j]

        cpe = etacoeff[p,j]
        origolinepoint = sqrt(cpe[1] / cpe[3])      # origoline touching a + bx + cx2 has x=sqrt(a/c) 
        origoline = Polynomial([0, cpe[2] + 2*cpe[3]*origolinepoint]) 

        segmentpoints = range(root, stop=maxturbdischarge[p,j], length=length(SEGMENT))
        origosegmentpoints = range(origolinepoint, stop=maxturbdischarge[p,j], length=length(SEGMENT))

        k_firstseg[p,j] = 0.3

        for s in SEGMENT
            x, ox = segmentpoints[s], origosegmentpoints[s]
            k_segcoeff[p,j,s] = dPE(x)
            m_segcoeff[p,j,s] = Eff_discharge(x) - k_segcoeff[p,j,s] * x
            k_segcoeff_origo[p,j,s] = dPE(ox)
            m_segcoeff_origo[p,j,s] = (s == 1) ? 0.0 : Eff_discharge(ox) - k_segcoeff_origo[p,j,s] * ox
        end

        end_rampseg[p,j] = m_segcoeff[p,j,1] / (k_firstseg[p,j] - k_segcoeff[p,j,1])
        end_rampseg[p,j] = (end_rampseg[p,j] > maxturbdischarge[p,j]) ? 0.0 : end_rampseg[p,j]
        end_zeroseg[p,j] = -m_segcoeff[p,j,1] / k_segcoeff[p,j,1]
        end_zeroseg_poly[p,j] = root
        end_origoseg[p,j] = origolinepoint

        bigM[p,j] = maxturbdischarge[p,j] - min(0, m_segcoeff[p,j,1])
        
        if !silent && i <= 0
            plotly()
            d = 0:.01:maxturbdischarge[p,j]

            plot(d, max.(k_firstseg[p,j] * d, Eff_discharge.(d)), lw=3)
            plot!(d, origoline.(d), lw=3)
            scatter!(p_discharge, p_pe, xlim=(0,Inf), ylim=(0,Inf), legend=nothing, xlabel="Discharge", ylabel="Discharge * Efficiency", title="Effective Discharge", tickfont=14, legendfont=14, guidefont=16, titlefont=18, size = (700, 500)) |> display
            #scatter!(dc, dc.*eta, label="Data", markersize=0.4, alpha=0.4) |> display

            y = [(dd <= end_rampseg[p,j]) ? k_firstseg[p,j] * dd : minimum(k_segcoeff[p,j,s] * dd + m_segcoeff[p,j,s] for s in SEGMENT) for dd in d]
            # y = [minimum(k_segcoeff_origo[p,j,s] * dd + m_segcoeff_origo[p,j,s] for s in SEGMENT) for dd in d]
            plot(d, y, lw=3)
            plot!(d, origoline.(d), lw=3)
            for s in SEGMENT
                plot!(d, k_segcoeff[p,j,s] .* d .+ m_segcoeff[p,j,s], line=(1, :lightgray))
                 plot!(d, k_segcoeff_origo[p,j,s] .* d .+ m_segcoeff_origo[p,j,s], line=(1, :lightgray))
            end
            scatter!(p_discharge, p_pe, xlim=(0,Inf), ylim=(0,Inf), legend=nothing, xlabel="Discharge", ylabel="Discharge * Efficiency", title="Effective Discharge", tickfont=14, legendfont=14, guidefont=16, titlefont=18, size = (700, 500)) |> display

            plot(d, max.(k_firstseg[p,j], Eff_discharge.(d)./d), lw=3)
            plot!(d, origoline.(d)./d, lw=3)
            #for s in SEGMENT
            #    plot!(d, k_segcoeff[p,j,s] .+ m_segcoeff[p,j,s] ./ d, line=(1, :lightgray))
            #    plot!(d, k_segcoeff_origo[p,j,s] .* d .+ m_segcoeff_origo[p,j,s], line=(1, :lightgray))
            #end
            scatter!(p_discharge, p_eta, xlim=(0,p_discharge[end].*1.05), ylim=(0,1), legend=nothing, xlabel="Discharge", ylabel="Efficiency", title="Efficiency", tickfont=14, legendfont=14, guidefont=16, titlefont=18, size = (700, 500)) |> display
            #scatter!(dc, eta, label="Data", markersize=0.4, alpha=0.4, xlim=(0,p_discharge[end].*1.05), ylim=(0,1), legend=:outertopright, size=(1000,600), title="$p $j eta") |> display
        end

        if !silent && i <= 1
            plotly()
            d = 0:.01:maxturbdischarge[p,j]
            dd = 0:.01:minimum(roots(Eff_discharge))

            # plotting effective discharge (discharge*eta)
            a = plot()
            for s in SEGMENT
                plot!(a, d, k_segcoeff[p,j,s] .* d .+ m_segcoeff[p,j,s], line=(2, :dodgerblue2, :dot), linealpha=1, label=(s<=1 ? "model A:MIP" : nothing))
            end
            for s in SEGMENT
                plot!(a, d, k_segcoeff_origo[p,j,s] .* d .+ m_segcoeff_origo[p,j,s], line=(2, :gold, :solid), la=1, label=(s<=1 ? "model B:L and C" : nothing))
            end
            plot!(a, d, origoline.(d), line=(7, :darkorange1, :dash), label="model D")
            plot!(a, d, Eff_discharge.(d), line=(7, :grey30), label="Typical turbine")
            plot!(a, dd, zeros(length(dd)), line=(7, :grey30, :solid), label=nothing)
            plot!(a, d, (d.<end_origoseg[p,j]).*k_segcoeff_origo[p,j,1].*d + (d.>=end_origoseg[p,j]).*Eff_discharge.(d), line=(7, :yellowgreen, :dot), label="model B")
            plot!(a, d, max.(k_firstseg[p,j] * d, Eff_discharge.(d)), line=(7, :dodgerblue2, :dash), label="model A")
            
            plot!(a, xlim=(0,Inf), ylim=(0,85), legend=:bottomright, xlabel="Discharge [m3/s]", ylabel="Effective discharge [m3/s]",
             tickfont=16, legendfont=16, guidefont=18, titlefont=20, size = (800, 500), automargin=true, left_margin=-5mm) |> display # title="Effective Discharge",

            # plotting efficiency
            ps=1
            b = plot()
            for s in SEGMENT
                plot!(b, d, k_segcoeff[p,j,s] .+ m_segcoeff[p,j,s] ./ d, line=(2, :dodgerblue2, :dot), linealpha=1, label=nothing) #(s<=1 ? "model A:MIP" : nothing))
            end 
            for s in SEGMENT
                plot!(b, d, k_segcoeff_origo[p,j,s] .+ m_segcoeff_origo[p,j,s]./d, line=(2*ps, :gold, :solid), la=1, label=nothing) #(s<=1 ? "B:L and C" : nothing))
            end
            plot!(b, d, origoline.(d)./d, line=(7*ps, :darkorange1, :dash), label=nothing) #"model D")
            plot!(b, d, Eff_discharge.(d)./d, line=(7*ps, :grey30, :solid), label=nothing) #"Typical turbine")
            plot!(b, dd, zeros(length(dd)), line=(7*ps, :grey30, :solid), label=nothing)
            plot!(b, d, max.(k_firstseg[p,j], Eff_discharge.(d)./d), line=(7*ps, :dodgerblue2, :dash), label=nothing) #"model A")
            plot!(b, d, (d.<end_origoseg[p,j]).*k_segcoeff_origo[p,j,1] + (d.>=end_origoseg[p,j]).*Eff_discharge.(d)./d, line=(7*ps, :yellowgreen, :dot), label=nothing) #"model B")
            
            
            
            
            plot!(b, xlim=(0, p_discharge[end].*1.05), ylim=(0,1), legend=nothing, xlabel=nothing, ylabel="Efficiency",
            tickfont=16*ps, legendfont=16*ps, guidefont=18*ps, titlefont=20*ps, size = (800*ps, 500*ps), automargin=true, left_margin=10mm) |> display #title="Efficiency",


            plot(b, a, layout=(2), size=(1200, 500), legend=:bottomright) |> display#, link=:x) |> display #title="Efficiency",

            #println("Maxdischarge (plant): ", sum(maxplantdischarge))
            #println("Maxdischarge (turb): ", sum(values(maxturbdischarge)))

        end
    end 

     #plotting approximation error
     #= for (i, p) in enumerate(PLANT), j in TURBINE[p]
        plants[i].capacity == 0 && continue
        d = 0:.01:maxturbdischarge[p,j]

        p_discharge = [etapoints[p,j][i].d for i in 1:length(etapoints[p,j])]
        p_eta = [etapoints[p,j][i].e for i in 1:length(etapoints[p,j])]
        p_pe = p_discharge .* p_eta
        Eff_discharge = Polynomials.fit(p_discharge, p_pe, 2)  # fit 2nd degree polynomial to Eff_discharge points (not eta!)
    
        meaneffectivedischarge = meanturbeta[p,j]*meanturbdischarge[p,j]
        maxeffectivedischarge = findmax(Eff_discharge.(d))[1]
        #println(" maxeffectivedischarge= ",  maxeffectivedischarge)
        E_range = 0:maxeffectivedischarge
        H_range = minhead[p]:maxhead[p]
        H_diff = H_range[end]-H_range[1]
        er = findmax([abs(E*H-E*meanhead[p]+meaneffectivedischarge*H-meanhead[p]*meaneffectivedischarge)*dens*grav*WtoMW for E in E_range, H in H_range])
        er2 = findmax([abs(E*H-E*maxhead[p])*dens*grav*WtoMW for E in E_range, H in H_range])
        Terror[p,j] = er[1]/(E_range[er[2][1]]*H_range[er[2][2]]*dens*grav*WtoMW)*100
        Cerror[p,j] = er2[1]/(E_range[er2[2][1]]*H_range[er2[2][2]]*dens*grav*WtoMW)*100
        #println("Plant: $p, Turbine: $j")
        #println("Maximum effective discharge= ", E_range[end])
        #println("Head variation= ", H_diff)
        #println("Maximum error Taylor= ", Terror[p,j])
        #println("Maximum error constant head= ", Cerror[p,j])

        #=     taylorerror = [i for i in values(error)]
        println(taylorerror)
        Cerror = [i for i in values(error2)]
        turbines = ["$p, $j" for (p,j) in keys(error2)]
        models = ["Taylor", "Constant head"]
        cats = repeat(turbines, inner=length(models))
        println(cats)
        groups = repeat(models, outer=length(turbines))
        println(groups)
        groupedbar(cats, [taylorerror Cerror; [1 1]]; group=groups) |> display =#
    end
    plotly()
    a = bar(collect(values(Terror)), label="Taylor", xlabel="Turbines", ylabel="Maximum error [% of prod that hour]")
    b = bar(collect(values(Cerror)), label="Constant head", fillalpha=0.5, xlabel="Turbines")  
plot(a,b,layout=2, link=:y, size=(1000, 500)) |> display  =#



    return etacoeff, maxturbdischarge, maxturbeta, meanturbdischarge, meanturbeta, 
            AxisArray(maxplantdischarge, PLANT), k_firstseg, end_rampseg, end_zeroseg, end_zeroseg_poly, end_origoseg,
            k_segcoeff, m_segcoeff, k_segcoeff_origo, m_segcoeff_origo, bigM 
end
