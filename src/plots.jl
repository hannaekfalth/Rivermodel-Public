using Plots.PlotMeasures

export plot_production, plot_production_plant, plot_head_plant, plot_spill_plant,
        plot_flow_plant, plot_inflow, plot_inflow_plant

const opt = (tickfont=18, legendfont=18, guidefont=20, titlefont=22, 
            legend = :outertopright, size = (1000, 500), automargin=true, left_margin=10mm)

const metricsopt = (xlabel="Computational time [s]", xaxis=:log, 
                xticks=(10 .^ range( 0, 6, length = 7), 10 .^ range( 0, 6, length = 7)), 
                minorgrid=true, minorticks=1000, top_margin=5mm)
        
function plot_production(years, runs, params, layouts)
    twoweeks = 876:1212
    plotly()
    for y in years, run in runs
        @unpack name, year, TIME, power = run
        year != y && continue
        letter, color, markershape, linestyle, _ = layouts[name]
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        power_time_series = sumdrop(power, dims=2)
        hist_prod_series = sumdrop(params[y].hist_prod, dims=2)[TIME]
        a = plot(TIME, power_time_series;
            label = nothing, margin=10mm, titleloc = :right, #title = "Power Production $year",
            xlabel = "Hour",  yticks=0:200:1000, ylabel = "MWh/h", linestyle, color, opt...)
        b = plot(TIME[twoweeks], power_time_series[twoweeks];
            label = letter, linestyle, linewidth = 2, xticks = 900:100:1300,
            xlabel = "Hour", color, opt..., left_margin=0mm, yformatter=_->"")
        plot!(a, TIME, hist_prod_series; color = hist_color, linestyle = hist_linestyle, label = nothing)
        plot!(b, TIME[twoweeks], hist_prod_series[twoweeks], linestyle = hist_linestyle, linewidth = 2,
            color = hist_color, label = hist_letter)
        plot(a, b, layout=2, link=:both, legend = :outertopright, size = (1200, 500), ylim=(-50,1200)) |> display
    end
end

function plot_production_comparetoA(years, runs, params, layouts)
    twoweeks=876:1212
    plotly()
    linew=2*2
    for y in years, run in runs
        @unpack name, year, TIME, power = run
        year != y && continue
        indexA = findfirst(r -> r.year == y && contains(layouts[r.name][1], "A"), runs)
        indexA === nothing && continue
        letter, color, markershape, linestyle, _ = layouts[name]
        letterA, colorA, markershapeA, linestyleA, _ = layouts[runs[indexA].name]
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        power_time_series = sumdrop(power, dims=2)
        hist_prod_series = sumdrop(params[y].hist_prod, dims=2)[TIME]
        A_power_time_series = sum(runs[indexA].power, dims=2)

        a = plot(TIME, A_power_time_series;
            label = nothing, #title = "Power Production $y", titleloc = :right,
            xlabel = "Hour",  yticks=0:200:1200, ylabel = "MWh/h", linestyle=linestyleA, linewidth=linew/2,  color=colorA, opt...)
        b = plot(TIME[twoweeks], A_power_time_series[twoweeks];
            label = letterA, linestyle=linestyleA, linewidth = linew, xticks = 900:100:1300,
            xlabel = "Hour", color=colorA, left_margin=0mm, yformatter=_->"", opt...) 
        c = plot(TIME, power_time_series;
                label = nothing, top_margin=-5mm,
                xlabel = "Hour",  yticks=0:200:1200, ylabel = "MWh/h", linestyle, linewidth=linew/2, color, opt...)
        d = plot(TIME[twoweeks], power_time_series[twoweeks];
                label = letter, linestyle, linewidth = linew, xticks = 900:100:1300,
                xlabel = "Hour", color, left_margin=0mm, top_margin=-5mm, yformatter=_->"", opt...)
        plot!(a, hist_prod_series, color=hist_color, linestyle=hist_linestyle, linewidth=linew/2, label = nothing)
        plot!(b, TIME[twoweeks], hist_prod_series[twoweeks], linestyle=hist_linestyle, linewidth = linew,
                color=hist_color, label = nothing)
        plot!(c, hist_prod_series, color=hist_color, linestyle=hist_linestyle, linewidth=linew/2, label = nothing)
        plot!(d, TIME[twoweeks], hist_prod_series[twoweeks], linestyle=hist_linestyle, linewidth = linew,
            color=hist_color, label = hist_letter)
        plot(a, b, c, d, layout=4, link=:both, ylim=(-50,1100), size = (1000*2, 500*2)) |> display
    end
end

function taylor_eval(years, runs, params, layouts)
    plotly()
    comparingwithtaylor = ["A", "B"]
    twoweeks=876:1212
    RMSD = @aa 0.0, comparingwithtaylor, Symbol.(years)
    b=2
    for y in years
        a = plot()
        for m in ["A", "B:L"]#, "B:L", "C", "D", "E"]
            runindex = findfirst(r -> r.year == y && layouts[r.name][1] == m, runs)
            runindex === nothing && continue
            power_time_series = sum(runs[runindex].power, dims=2)
            letter, color, markershape, linestyle, _ = layouts[runs[runindex].name]
            plot!(a, params[y].TIME[twoweeks], power_time_series[twoweeks]; xlabel = "Hour", ylabel = "MWh/h", 
                label = m, linestyle, color, linewidth = 3*b, opt..., tickfont=14*b, legendfont=16*b, guidefont=16*b, titlefont=16*b, titleloc=:center, 
                legend = :outertopright, size = (700*b, 400*b), automargin=true, left_margin=10mm,)
        end
        display(a)
        for m in comparingwithtaylor
            runindex = findfirst(r -> r.year == y && layouts[r.name][1] == m, runs)
            runindex_BT = findfirst(r -> r.year == y && layouts[r.name][1] == "B:L", runs)
            (runindex === nothing || runindex_BT === nothing) && continue
            power_time_series = sum(runs[runindex].power, dims=2)
            T_power_time_series = sum(runs[runindex_BT].power, dims=2)
            RMSD[m,Symbol(y)] = rmsd(power_time_series, T_power_time_series; normalize=false)/maximum(power_time_series)*100
        end
    end
    legends = ["Comparing B:L to $m" for m in comparingwithtaylor]
    groupedruns = repeat(legends, outer=length(years)) |> CategoricalArray
    namedcolors = [val[2] for (key,val) in layouts if val[1] in comparingwithtaylor]
    fillalphas = [val[5] for (key,val) in layouts if val[1] in comparingwithtaylor]
    groupedbar(repeat(string.(years), inner=length(comparingwithtaylor)), RMSD;
        group = groupedruns, xlabel = "Year", ylabel = "RMSD on hourly prod. [%]", # % of max. prod. of A resp. B  
        left_margin=10mm, opt...,
        color = repeat(namedcolors, outer=length(years)),
        linecolor = :match,
        fillalpha = repeat(fillalphas, outer=length(years)),
    ) |> display
end

function plot_production_plant(years, runs, params, layouts; runletters, plants)
    plotly()
    for y in years, run in runs
        @unpack name, year, TIME, power = run
        year != y && continue
        contains(name, "aggregated") && continue
        letter, color, markershape, linestyle, _ = layouts[name]
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        !isempty(runletters) && !(letter in runletters) && continue
        for p in params[y].PLANT
            !isempty(plants) && !(p in plants) && continue
            plot(TIME, vec(power[:,p]); label = letter, color)
            plot!(TIME, params[y].hist_prod[:,p]; label = hist_letter, color=hist_color,
                title = "Power Production $y - $p", xlabel = "Hour", ylabel = "MWh/h") |> display
        end
    end 
end

function plot_production_turbine(years, runs, params, layouts; runletters, plants)
    plotly()
    for y in years, run in runs
        @unpack name, year, TIME, power = run
        year != y && continue
        contains(name, "aggregated") && continue
        letter = layouts[name][1]
        !isempty(runletters) && !(letter in runletters) && continue
        for p in params[y].PLANT
            !isempty(plants) && !(p in plants) && continue
            plot()
            for j in params[y].TURBINE[p]
                plot!([run.power_turbines[t,p,j] for t in TIME]; label = "Turbine $j")
            end
            plot!(; title = "Power Production $letter - $y - $p", xlabel = "Hour", ylabel = "MWh/h", opt...) |> display
        end
    end 
end

function plot_power_duration(years, runs, params, layouts)
    plotly()
    for y in years
        plot()
        
        for run in runs
            @unpack name, year, power = run
            year != y && continue
            sortedpower = sort(sumdrop(power, dims=2), rev=true)
            letter, color, markershape, linestyle, _ = layouts[name]
            plot!(sortedpower; linewidth = 5, label=letter, color, linestyle)
        end
        sortedhistpower = sort(sumdrop(params[y].hist_prod, dims=2), rev=true)
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        plot!(sortedhistpower; title = "PDC - $y", xlabel = "Hour", ylabel = "MWh/h", 
                linewidth=3, label=hist_letter, color=hist_color, linestyle=hist_linestyle, opt...) |> display

       
        plot()
        plants = getproperty.(params[y].plants, :capacity) .> 0
        plants[params[y].PLANT .== :Bergnäs] .= false    # avoid 400m data error
        for run in runs
            @unpack name, year, head = run
            year != y && continue
            contains(name, "aggregated") && continue
            sortedhead = sort(meandrop(head[:,plants], dims=2), rev=true)
            letter, color, markershape, linestyle, _ = layouts[name]
            plot!(sortedhead; linewidth=3, label=letter, color, linestyle)
        end
        
        sortedhisthead = sort(meandrop(params[y].hist_head[:,plants], dims=2), rev=true)
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        plot!(sortedhisthead; title = "HDC - $y", xlabel = "Hour", ylabel = "m", 
                linewidth=3, label=hist_letter, color=hist_color, linestyle=hist_linestyle, opt...) #|> display
    end
end

function plot_power_duration_plant(years, runs, params, layouts; plants)
    plotly()
    for y in years, (i, p) in enumerate(params[y].PLANT)
        !isempty(plants) && !(p in plants) && continue
        plot()
        for run in runs
            @unpack name, year, power = run
            year != y && continue
            contains(name, "aggregated") && continue
            sortedpower = sort(power[:,i], rev=true)
            letter, color, markershape, linestyle, _ = layouts[name]
            plot!(sortedpower; linewidth = 3, label=letter, color, linestyle)
        end
        sortedhistpower = sort(params[y].hist_prod[:,i], rev=true)
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        plot!(sortedhistpower; title = "PDC $p - $y", xlabel = "Hour", ylabel = "MWh/h", 
                linewidth=3, label=hist_letter, color=hist_color, linestyle=hist_linestyle, opt...) |> display

        plot()
        for run in runs
            @unpack name, year, head = run
            year != y && continue
            contains(name, "aggregated") && continue
            sortedhead = sort(head[:,i], rev=true)
            letter, color, markershape, linestyle, _ = layouts[name]
            plot!(sortedhead; linewidth=3, label=letter, color, linestyle)
        end
        sortedhisthead = sort(params[y].hist_head[:,i], rev=true)
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        plot!(sortedhisthead; title = "HDC $p - $y", xlabel = "Hour", ylabel = "m", 
                linewidth=3, label=hist_letter, color=hist_color, linestyle=hist_linestyle, opt...) |> display
    end
end

function plot_prod_price(years, runs, params, layouts)
    plotly()
    for y in years
        price = params[y].spot_price
        plot()
        for run in runs
            @unpack name, year, power = run
            year != y && continue
            letter, color, markershape, linestyle, _ = layouts[name]
            scatter!(price, sumdrop(run.power, dims=2); label=letter, 
                markersize=3, markerstrokewidth=0, alpha=0.1, color)
        end
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        scatter!(price, sumdrop(params[y].hist_prod, dims=2); label=hist_letter,
                color=hist_color, markersize=3, markerstrokewidth=0, alpha=0.1,
                title = "Production vs price - $y", xlabel = "Price [SEK/MWh]", ylabel = "Production [MWh/h]",
                opt...) |> display
    end
end

function plot_prod_price_article(years, runs, params, layouts)
    plotly()
    A_name = [k for (k,v) in layouts if v[1] == "A"][1]
    E_name = [k for (k,v) in layouts if v[1] == "E"][1]
    HD_name = [k for (k,v) in layouts if v[1] == "Historical Data"][1]
    a = plot()
    b = plot()
    for y in years
        price = params[y].spot_price
        plot()
        for run in runs
            @unpack name, year, power = run
            year != y && continue
            if name == A_name || name == E_name
            letter, color, markershape, linestyle, _ = layouts[name]
            scatter!(b, price, sumdrop(run.power, dims=2); label=letter, 
                markersize=2*2, markerstrokewidth=0, alpha=0.5, color, xlabel = "Price [SEK/MWh]",  opt..., leftmargin=-5mm)
            end
        end
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        scatter!(a, price, sumdrop(params[y].hist_prod, dims=2); label=hist_letter,
                color=hist_color, markersize=2*2, markerstrokewidth=0, alpha=0.5,
                #title = "Production vs price - $y",  
                ylabel = "Production [MWh/h]", xlabel = "Price [SEK/MWh]",
                opt..., rightmargin=-5mm ) 
        plot(a, b, layout=2, link=:y, size=(800*2,400*2), legend=(0.8,0.4)) |> display
    end
end

function plot_prod_price_plant(years, runs, params, layouts; plants)
    plotly()
    for y in years, (i, p) in enumerate(params[y].PLANT)
        !isempty(plants) && !(p in plants) && continue
        price = params[y].spot_price
        plot()
        for run in runs
            @unpack name, year, power = run
            year != y && continue
            contains(name, "aggregated") && continue
            letter, color, markershape, linestyle, _ = layouts[name]
            scatter!(price, power[:,i]; label=letter, 
                markersize=1.5, markerstrokewidth=0, alpha=0.1, color)
        end
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        scatter!(price, params[y].hist_prod[:,i]; label=hist_letter,
                color=hist_color, markersize=3, markerstrokewidth=0, alpha=0.1,
                title = "Production vs price - $p - $y", xlabel = "Price [SEK/MWh]", ylabel = "Production [MWh/h]",
                opt...) |> display
    end
end

function plot_prod_price_all(years, runs, params, layouts)
    plotly()
    for y in years
        price = params[y].spot_price
        aa=0.3  #alpha
        ms=2    #markersize
        msw=0   #markerstrokewidth

        a,b,c,d,e,f,g,h = (plot() for _ =1:8)
        plots = [a,b,c,d,e,f,g,h]
        plot()
        i = 0
        for run in runs
            @unpack name, year, power = run
            year != y && continue
            i += 1
            totpower = sumdrop(power, dims=2)
            letter, color, markershape, linestyle, _ = layouts[name]
            scatter!(plots[i], price, totpower; label=letter, title=letter, 
                markersize=ms, markershape, markerstrokewidth=msw, alpha=aa, color)
        end
        totpower = sumdrop(params[y].hist_prod, dims=2)
        hist_letter, hist_color, hist_markershape, hist_linestyle, _ = layouts["Historical Data"]
        scatter!(plots[end], price, totpower; label=hist_letter, title=hist_letter,
                color=hist_color, markersize=ms, markerstrokewidth=msw, alpha=aa)
        
        plot(a,b,c,d,e,f,g,h, link = :both, layout=@layout([° °; ° °; ° °; ° °]), xlabel = "Price [SEK/MWh]", ylabel = "Production [MWh/h]", tickfont=12*2, legendfont=14*2, guidefont=14*2, titlefont=14*2, 
        legend = false, size = (1200*2, 1500*2)) |> display #empty instead of a letter to get empty place in figure
    end
end

function plot_plant_investigation(years, runs, params, layouts)
    plotly()
    for y in years
        @unpack hours, production, forebay_level, tail_level, discharge = params[y].vars

        year = Dates.year.(hours)
        yr = (year .>= 2015) .& (year .<= 2020)

        priceorder = sortperm(params[y].spot_price, rev=true)

        kvist = findfirst(params[y].PLANT .== :Kvistforsen)
        up = @view forebay_level[yr, kvist]
        quick_clean!(up, mean(up[up.>0]))
        low = @view tail_level[yr, kvist]
        quick_clean!(low, mean(low[low.>0]))

        nplants = length(params[y].PLANT)
        dcmax = zeros(nplants)

        for (i, p) in enumerate(params[y].PLANT)
            a = plot()
            b = plot()
            c = plot()
            d = plot()
            for (name, q) in runorder
                if contains(name, "aggregated")
                    nothing
                elseif q==y
                    plot!(a, rundata[name, y].power[:,i],    #plant production
                        label = name, color=colors[name], legend=false)
                    scatter!(b, rundata[name, y].power[:,i][priceorder]; label = name, color=colors[name], markersize=2, markerstrokewidth=0, alpha=0.5)
                    #plot!(b, sort!(rundata[name, y].power[:,i], rev=true); label = name, color=colors[name], linewidth = 2.2)
                end
                if contains(name, "ta-function") && q==y
                    params[y].plants[i].capacity == 0 && continue
        
                    head = forebay_level[yr, i] - tail_level[yr, i]
                    dc0 = discharge[yr, i]
                    calc_elec = 9.81*997*1e-6 * dc0 .* head
                    eta0 = production[yr, i] ./ calc_elec
                    f = (eta0 .> 0) .& (eta0 .< 1) .& (forebay_level[yr, i] .> 0) .& (tail_level[yr, i] .> 0) .&
                                (dc0 .> 0) .& (production[yr, i] .> 0)
                    dc, eta = dc0[f], eta0[f] # discharge and eta per hour for plant i
                    dcmax[i] = quantile(dc, 0.999)
        
                    
        
                    c=scatter(dc, eta, xlim=[0, dcmax[i]*1.1], ylim=[0,1], label="Data", markersize=0.8, markerstrokewidth=0, alpha=0.2)
        
                    for j in params[y].TURBINE[p]
                        #plot!(c, [rundata["Eta-function_NLP", y].power_turbines[t,p,j] for t in 1:length(params[y].TIME)],
                        #label = "Turbine $j", color=turbcolors[j+6], title = "Production per turbine", xlabel = "Hour", ylabel = "MWh/h",
                        #tickfont=14, guidefont=16, titlefont=18)
                        

                        etafunction(x) = (params[y].etacoeff[p, j][1] + params[y].etacoeff[p, j][2]*x + params[y].etacoeff[p, j][3]*x^2)
                        plot!(c, etafunction, xlim=[0, params[y].maxdischarge[p]*1.1], ylim=[0,1], color=turbcolors[j+6], title="eta(discharge)", label="Turbine $j",
                        xlabel = "Discharge [m3/s]", ylabel = "efficiency", tickfont=12, guidefont=12, titlefont=14)
                        
                        plot!(d, [rundata["Eta-function_NLP", y].discharge_turbines[t,p,j] for t in 1:length(params[y].TIME)],
                        label = "Turbine $j", color=turbcolors[j+6], title = "Discharge per turbine", xlabel = "Hour", ylabel = "m3/s",
                        tickfont=12, guidefont=12, titlefont=14)
                    end
                else nothing
                end
            end
            plot!(a, params[y].TIME, params[y].hist_prod[:,p], label = "Historical Data", color=colors["Historical Data"],
                title = "Power Production", xlabel = "Hour", ylabel = "MWh/h",
                tickfont=12, guidefont=12, titlefont=14)
            #plot!(b, sort(params[y].hist_prod[:,p], rev=true), label = "Historical Data", color=colors["Historical Data"],
            #    linewidth = 2.2, title = "Power Duration Curve", xlabel = "Hour", ylabel = "MWh/h",
            #    tickfont=12, guidefont=12, titlefont=14)
            scatter!(b, params[y].hist_prod[:,p][priceorder], label = "Historical Data", color=colors["Historical Data"],
                markersize=2, markerstrokewidth=0, alpha=0.5, title = "Power Duration Curve", xlabel = "Hour", ylabel = "MWh/h",
                tickfont=12, guidefont=12, titlefont=14)
            title=plot(title = "$p - $y", grid = false, showaxis = false, titlefont=18, bottom_margin = -10mm)
            
            
            l = @layout [title{0.01h}; a b; d c]
            plot(title, a, b, d, c, layout=l, link = :xaxis, size = (1900, 900), legend = :outertopright, 
            margin = 10mm, left_margin=10mm, legendfont=13) |> display
        end
    end 
end

function metricsvsCT(years, runs, params, layouts)
    nrmetrics=8
    RMSD_CT, RMSDramp_CT, PC_CT, Rev_CT, Ramp_CT, Max_CT, Min_CT, meanerror = (plot() for _ =1:nrmetrics)
    plots = [RMSD_CT, RMSDramp_CT, PC_CT, Rev_CT, Ramp_CT, Max_CT, Min_CT, meanerror]
    b = 1 # figure scaling setting
    
    for (k, year) in enumerate(years)
        @unpack PLANT, TIME, plant, hist_prod, hist_head, spot_price = params[year]
        hist_prod_tot = sum(hist_prod[:, p] for p in PLANT)
        ramping_hist_series = abs.(diff(hist_prod_tot)) 
        indexA = findfirst(r -> r.year == year && contains(layouts[r.name][1], "A"), runs)
        indexA === nothing && continue
        Apower = sumdrop(runs[indexA].power, dims=2)
        A_ramping_series = abs.(diff(Apower)) 
        runs_year = [run for run in runs if run.year == year && run.name in keys(layouts)]
        RMSD, RMSDramp, PC, CT, Rev, Ramp, Maxprod, Minprod, me = (Vector{Float64}(undef, length(runs_year)) for _ = 1:nrmetrics+1)
        metrics = [RMSD, RMSDramp, PC, Rev, Ramp, Maxprod, Minprod, me]

        for (i, run) in enumerate(runs_year)
            totalpower = sumdrop(run.power, dims=2)
            ramping_series = abs.(diff(totalpower))

            RMSD[i] = rmsd(Apower, totalpower; normalize=false)#./maximum(hist_prod_tot).*100
            RMSDramp[i] = rmsd(A_ramping_series, ramping_series; normalize=false)#./maximum(ramping_hist_series).*100
            PC[i] = cor(hist_prod_tot, totalpower)
            Rev[i] = (run.profit/(sum(hist_prod_tot.*spot_price)/1e6)-1)*100
            Ramp[i] = ((mean(sort(ramping_series, rev=true)[1:round(Int, 0.01*length(ramping_series))])./mean(sort(ramping_hist_series, rev=true)[1:round(Int, 0.01*length(ramping_hist_series))]))-1)*100
            Maxprod[i] = (maximum(totalpower)/maximum(hist_prod_tot)-1)*100
            Minprod[i] = (minimum(totalpower)/minimum(hist_prod_tot)-1)*100
            me[i] = mean(abs.(Apower-totalpower))
            CT[i] = run.solvetime

            letter, color, markershape, linestyle, _ = layouts[run.name]

            for m in 1:nrmetrics
                scatter!(plots[m], [CT[i]], [metrics[m][i]]; color, markershape, markersize=6*b, markeralpha=0.8, label=(k>1 ? nothing : letter)) #label=(k>1 ? nothing : names[name])
            end
        end
        #for m in 1:nrmetrics 
        #plot!(plots[m], CT, metrics[m], linecolor=:black, label=(m>1 ? nothing : "Same year")) #label=(k>1 ? nothing : "Same year")
        #end
    end
    
    scatter!(RMSD_CT, #title="b)",
                    ylabel="RMSD on hourly production [MWh/h]", ylim=[0,600], xlabel="Computational time [s]", xaxis=:log, 
                    xticks=(10 .^ range( 0, 6, length = 7), 10 .^ range( 0, 6, length = 7)), 
                    tickfont=12*b, legendfont=14*b, guidefont=14*b, titlefont=14*b, titleloc=:center, 
                    legend = :outertopright, size = (500*b, 400*b), automargin=true, left_margin=10mm,
                    minorgrid=true, minorticks=1000, top_margin=5mm, gridlinewidth=0.1, gridalpha = 1, gridcolor=:grey80, thickness_scaling=4) #|> display
                    display(RMSD_CT)
    scatter!(RMSDramp_CT, #title="c)",
                    ylabel="RMSD on hourly ramping [MWh/h]", ylim=[0,300], xlabel="Computational time [s]")
    scatter!(PC_CT, #title="b)",
                    ylabel="Correlation with historical data", ylim=[-1,1])
    scatter!(Rev_CT, title="a)",
                    ylabel="Deviation in revenue [%]", ylim=[0,15])
    scatter!(Ramp_CT, title="d)",
                    ylabel="Deviation in maximum ramping [%]", ylim=[0,300])
    scatter!(Max_CT, title="e)",
                    ylabel="Deviation in maximum production [%]", ylim=[0,25])
    scatter!(Min_CT, title="f)",
                    ylabel="Deviation in minimum production [%]", ylim=[-110,110]) 
    scatter!(meanerror,
                    ylabel="Mean error [MWh/h]", ylim=[0,600], xlabel="Computational time [s]", xaxis=:log, 
                    xticks=(10 .^ range( 0, 6, length = 7), 10 .^ range( 0, 6, length = 7)), 
                    tickfont=12*4, legendfont=14*4, guidefont=14*4, titlefont=14*4, titleloc=:center, 
                    legend = :outertopright, size = (500*4, 400*4), automargin=true, left_margin=10mm,
                    minorgrid=true, minorticks=1000, top_margin=5mm, gridlinewidth=0.1, gridalpha = 1, gridcolor=:grey80) |> display

    #title=plot(title = "$p - $y", grid = false, showaxis = false, titlefont=18, bottom_margin = -10mm)        
    #=plot(Rev_CT, RMSD_CT, RMSDramp_CT, Ramp_CT, Max_CT, Min_CT, layout=@layout([° °; ° °; ° °]),
                    link=:x, #title = ["a","b","c","d","e","f"],
                    xlabel="Computational time [s]", xaxis=:log, 
                    xticks=(10 .^ range( 0, 6, length = 7), 10 .^ range( 0, 6, length = 7)), 
                    tickfont=12, legendfont=14, guidefont=14, titlefont=14, titleloc=:center, 
                    legend = :outertopright, size = (1500, 1500), automargin=true, left_margin=10mm,
                    minorgrid=true, minorticks=1000, top_margin=5mm) |> display 
    =#

    plot(RMSD_CT, RMSDramp_CT, layout=@layout([°; °]),
                    link=:x, #title = ["a","b","c","d","e","f"],
                    xlabel="Computational time [s]", xaxis=:log, 
                    xticks=(10 .^ range( 0, 6, length = 7), 10 .^ range( 0, 6, length = 7)), 
                    tickfont=12, legendfont=14, guidefont=14, titlefont=14, titleloc=:center, 
                    legend = :outertopright, size = (500, 700), automargin=true, left_margin=10mm,
                    minorgrid=true, minorticks=1000, top_margin=5mm, gridlinewidth=0.1, gridalpha = 1, gridcolor=:grey80)
end

function plotprofit(years, origruns, params, layouts)
    runs = copy(origruns)
    for (i, run) in enumerate(runs)
        if i < length(runs)
            letter, nextletter = layouts[run.name][1], layouts[runs[i+1].name][1]
            recalc = endswith(letter, "(r)") ? -1.0 : endswith(nextletter, "(r)") ? runs[i+1].profit : 0.0
        else
            recalc = 0.0
        end
        runs[i] = (run..., recalc=recalc)
    end
    filter!(r -> r.recalc >= 0, runs)

    plotly()
    runnames = unique([run.name for run in runs if run.name in keys(layouts)])
    runnames_hist = [runnames; "Historical Data"]
    letters_hist = [layouts[name][1] for name in runnames_hist]
    Rev = @aa 0.0, letters_hist, Symbol.(years)
    Rev_recalc = @aa 0.0, letters_hist, Symbol.(years)
    for run in runs
        letter = layouts[run.name][1]
        Rev[letter, Symbol(run.year)] = run.profit
        Rev_recalc[letter, Symbol(run.year)] = run.recalc
    end
    for y in years
        @unpack PLANT, TIME, hist_prod, spot_price = params[y]
        hist_prod_tot = sum(hist_prod[:, p] for p in PLANT)
        profit_hist = sum(hist_prod_tot[t]*spot_price[t] for t in TIME)/1e6
        Rev["Historical Data", Symbol(y)] = profit_hist
    end

    groupedruns = repeat(letters_hist, outer=length(years)) |> CategoricalArray # same hack again
    levels!(groupedruns, letters_hist)
    namedcolors = [layouts[name][2] for name in runnames_hist]
    fillalphas = [layouts[name][5] for name in runnames_hist]
    
    groupedbar(repeat(string.(years), inner=length(runnames_hist)), Rev; hover=Rev,
            group = groupedruns, xlabel = "Year", ylabel = "Revenue [Million SEK]", 
            #title = "Yearly revenue for the different models", 
            left_margin=10mm, opt...,
            size=(450+150*length(years),500),
            ylim=(0, 5+round(maximum(Rev), RoundUp, digits=-2)), #round(minimum(Rev), RoundDown, digits=-2)
            color = repeat(namedcolors, outer=length(years)),
            linecolor = :match,
            fillalpha = repeat(fillalphas, outer=length(years)),
            )

    barlines_flipped!(Rev_recalc) |> display
end

function articlemetrics(years, runs, params, layouts)
    plotly()
    runnames = unique([run.name for run in runs if run.name in keys(layouts)])
    runnames_hist = [runnames; "Historical Data"]
    letters = [layouts[name][1] for name in runnames]
    letters_hist = [layouts[name][1] for name in runnames_hist]
    Rev = @aa 0.0, letters_hist, Symbol.(years)
    Prod = @aa 0.0, letters_hist, Symbol.(years)
    RMSD = @aa 0.0, letters, Symbol.(years)
    RMSDplant = @aa 0.0, letters, Symbol.(years)
    RMSDramp = @aa 0.0, letters, Symbol.(years)

    warmstart = @aa NaN, letters_hist, Symbol.(years)   # warm starts given to NLP models
    ubound = @aa NaN, letters_hist, Symbol.(years)      # upper bounds from running MIP-model

    for year in years
        
        @unpack PLANT, TURBINE, TIME, plant, hist_prod, spot_price = params[year]

        hist_prod_tot = sum(hist_prod[:, p] for p in PLANT)
        ramping_hist_series = abs.(diff(hist_prod_tot)) 
        profit_hist = sum(hist_prod_tot[t]*spot_price[t] for t in TIME)/1e6
        
        for run in runs
            @unpack name, TIME, power = run
            run.year != year && continue
            letter, color, markershape, linestyle, _ = layouts[name]
            totalpower = sumdrop(run.power, dims=2)
            indexA = findfirst(r -> r.year == year && (layouts[r.name][1] == "A" || layouts[r.name][1] == "A_B:L" || layouts[r.name][1] == "A_C" || layouts[r.name][1] == "A_MIP"), runs)
            indexA === nothing && error("No A run for $year found in runs.")
            Apower = runs[indexA].power
            
            ramping_series = abs.(diff(totalpower)) 
            Rev[letter, Symbol(year)] = run.profit
            Prod[letter, Symbol(year)] = sum(totalpower) / 1e6
            RMSD[letter, Symbol(year)] = rmsd(sumdrop(Apower, dims=2), totalpower; normalize=false)/maximum(sumdrop(Apower, dims=2))*100
            if !contains(name, "aggregated") 
                rmsdplants = [rmsd(Apower[:,j], run.power[:,j]; normalize=true) for (j,p) in enumerate(PLANT) if plant[p].capacity > 0]
                RMSDplant[letter, Symbol(year)] = mean(rmsdplants)
            end
            RMSDramp[letter, Symbol(year)] = rmsd(ramping_hist_series, ramping_series; normalize=false)
            if run.warmstart !== nothing
                warmstart[letter, Symbol(year)] = run.warmstart
            end
            ubound[letter, Symbol(year)] = MIP_upper_bounds(letter, year)
        end

        Rev[layouts["Historical Data"][1], Symbol(year)] = profit_hist
        Prod[layouts["Historical Data"][1], Symbol(year)] = sum(hist_prod_tot) / 1e6
        
    end

    groupedruns = repeat(letters_hist, outer=length(years)) |> CategoricalArray # same hack again
    levels!(groupedruns, letters_hist)
    namedcolors = [layouts[name][2] for name in runnames_hist]
    fillalphas = [layouts[name][5] for name in runnames_hist]
    
    groupedbar(repeat(string.(years), inner=length(runnames_hist)), Rev; hover=Rev,
            group = groupedruns, xlabel = "Year", ylabel = "Revenue [Million SEK]", 
            #title = "Yearly revenue for the different models", 
            left_margin=10mm, opt...,
            size=(450*2+50*2*length(years),500*2),
            ylim=(0, 5+round(maximum([Rev; ubound]), RoundUp, digits=-2)), #round(minimum(Rev), RoundDown, digits=-2)
            color = repeat(namedcolors, outer=length(years)),
            linecolor = :match,
            fillalpha = repeat(fillalphas, outer=length(years)),
            )
    errorbars_flipped!(Rev, warmstart, ubound; marker_high=(:hline, :black, 6-length(years)/2)) |> display

    groupedbar(repeat(string.(years), inner=length(runnames_hist)), Prod; hover=Prod,
            group = groupedruns, xlabel = "Year", ylabel = "Annual production [TWh/year]", 
            #title = "Yearly revenue for the different models", 
            left_margin=10mm, opt...,
            size=(450*2+50*2*length(years),500*2),
            ylim=(0, round(maximum(Prod), RoundUp, digits=1) + 0.1), #round(minimum(Prod), RoundDown, digits=1) - 0.1
            color = repeat(namedcolors, outer=length(years)),
            linecolor = :match,
            fillalpha = repeat(fillalphas, outer=length(years)),
            ) |> display

    groupedruns2 = repeat(letters, outer=length(years)) |> CategoricalArray # same hack again
    levels!(groupedruns2, letters)
    namedcolors2 = [layouts[name][2] for name in runnames]
    fillalphas2 = [layouts[name][5] for name in runnames]

    groupedbar(repeat(string.(years), inner=length(runnames)), RMSDplant; hover=RMSDplant,
            group = groupedruns2, xlabel = "Year", ylabel = "RMSD on hourly plant production [%]", 
            #title = "Yearly revenue for the different models", 
            left_margin=10mm, opt...,
            size=(450*2+50*2*length(years),500*2),
            color = repeat(namedcolors2, outer=length(years)),
            linecolor = :match,
            ylim=[0,1],
            fillalpha = repeat(fillalphas2, outer=length(years)),) |> display

    groupedbar(repeat(string.(years), inner=length(runnames)), RMSD; hover=RMSD,
            group = groupedruns2, xlabel = "Year", ylabel = "RMSD on hourly prod. [% of max. prod. in A]", 
            #title = "Yearly revenue for the different models", 
            left_margin=10mm, opt...,
            size=(450*2+50*2*length(years),500*2),
            color = repeat(namedcolors2, outer=length(years)),
            linecolor = :match,
            ylim=[0,50],
            fillalpha = repeat(fillalphas2, outer=length(years)),) |> display
            
    groupedbar(repeat(string.(years), inner=length(runnames)), RMSDramp; hover=RMSDramp,
            group = groupedruns2, xlabel = "Year", ylabel = "RMSD on hourly ramping [MWh/h]", 
            #title = "Yearly revenue for the different models", 
            left_margin=10mm, opt...,
            size=(450*2+50*2*length(years),500*2),
            color = repeat(namedcolors2, outer=length(years)),
            linecolor = :match,
            ylim=[0,600],
            fillalpha = repeat(fillalphas2, outer=length(years)),) |> display

end  

function multiple_run_metrics(years, runs, params, layouts)
    plotly()
    allrunnames = unique([run.name for run in runs])
    RMSE = @aa 0.0, allrunnames, Symbol.(years)

    for year in years
        
        @unpack PLANT, TIME, plant, hist_prod, hist_head, spot_price = params[year]

        installed_capacity = sum(plant[p].capacity for p in PLANT)
        hist_prod_tot = sum(hist_prod[:, p] for p in PLANT)
        hist_head_tot = sum(hist_head[:, p] for p in PLANT if plant[p].capacity > 0)
        ramping_hist = maximum(abs.(diff(hist_prod_tot)))
        avg_ramping_hist = mean(abs.(diff(hist_prod_tot)))
        annual_prod_hist = sum(hist_prod)/1000000 #TWh
        profit_hist = sum(hist_prod_tot[t]*spot_price[t] for t in TIME)/1e6

        yearruns = [run for run in runs if run.year == year]

        profit, annual, ramping, avg_ramping, maxprod, minprod, highnr, lownr,
                RMSD, RMSD_head = (Vector{Float64}(undef, length(yearruns)) for _ = 1:10)
        b = Vector{Tuple{Float64, Float64, Float64}}(undef, length(yearruns))
        
        for (i, run) in enumerate(yearruns)
            totalpower = sumdrop(run.power, dims=2)
            profit[i] = run.profit
            #b[i] = bounds[layouts[name][1], Symbol(year)] 
            annual[i] = sum(totalpower)/1000000   #TWh
            ramping[i] = maximum(abs.(diff(totalpower)))
            avg_ramping[i] = mean(abs.(diff(totalpower)))
            maxprod[i] = maximum(totalpower)
            minprod[i] = minimum(totalpower)
            highnr[i] = count(x -> x >= 0.85*installed_capacity, totalpower)
            lownr[i] = count(x -> x <= 0.15*installed_capacity, totalpower)
            RMSD[i] = rmsd(hist_prod_tot, totalpower; normalize=false)
            RMSE[run.name, Symbol(year)] = RMSD[i]

            if contains(run.name, "aggregated")
                RMSD_head[i] = NaN
            else
                plants = [j for (j, p) in enumerate(PLANT) if plant[p].capacity > 0]
                #totalhead = sumdrop(rundata[name, year].head[:, plants], dims=2)
                #RMSD_head[i] = rmsd(hist_head_tot, totalhead; normalize=false)
            end
        end

        runnames = [[run.name for run in yearruns]; "Historical Data"]
        letters = [layouts[name][1] for name in runnames]
        namedcolors = [layouts[name][2] for name in runnames]
        fillalphas = [layouts[name][5] for name in runnames]

        metricnames = ["Annual Prod.", "Max Prod.", "Min Prod.", "Ramping Max.", "Ramping Avg."] #,"Profit"] 
        results = [ [annual./annual_prod_hist maxprod./maximum(hist_prod_tot) minprod./minimum(hist_prod_tot) ramping./ramping_hist avg_ramping./avg_ramping_hist];# profit./profit_hist];
                    [1 1 1 1 1]]
        # super ugly hack to fix groupedbar x-axis order bug
        # https://github.com/JuliaPlots/StatsPlots.jl/issues/291
        groupedruns = repeat(letters, outer=length(metricnames)) |> CategoricalArray
        levels!(groupedruns, letters)
        block = repeat(metricnames, inner=length(runnames))

        plotly()
        groupedbar(block, results*100; group=groupedruns,
            xlabel= "Metric", ylabel = "% of historical data", ylim=(0,280), 
            title = "Different models compared to historical data - $(year)", opt...,
            fillalpha = repeat(fillalphas, outer=length(metricnames)),
            color = repeat(namedcolors, outer=length(metricnames))) |> display

        println("\n\n\n\nYear: $year")
        println("Run order: ", letters[1:end-1])
        println()
        println("Revenue (Million SEK): ", round.(profit, digits=0)), 
        println("Bounds (Million SEK): ", b)
        println("Revenue (Million SEK - historical data: ", round.(profit_hist, digits=0)),
        println()
        println("Total power production - Model (TWh): ", round.(annual, digits=2))
        println("Total power production - Historical Data (TWh): ", round.(annual_prod_hist, digits=2))
        println()
        println("Weighted revenue - Model (SEK/MWh): ", round.(profit./annual, digits=0))
        println("Weighted revenue - Historical Data (SEK/MWh): ", round.(profit_hist./annual_prod_hist, digits=0))
        println()
        println("Deviation in production from historical data (%): ", round.((annual.-annual_prod_hist)/annual_prod_hist*100, digits=2))
        println("Deviation in revenue from historical data (%): ", round.((profit.-profit_hist)/profit_hist*100, digits=2))
        println("Deviation in weighted revenue from historical data (%): ", round.(((profit./annual).-(profit_hist./annual_prod_hist))/(profit_hist./annual_prod_hist)*100, digits=2))
        println()
        println("Maximum Ramping - Model (MWh/h): ", round.(ramping, digits=1))
        println("Maximum Ramping - Historical Data (MWh/h): ", round(ramping_hist, digits=1))
        println()
        println("Average Ramping - Model (MWh/h): ", round.(avg_ramping, digits=1))
        println("Average Ramping - Historical Data (MWh/h): ", round(avg_ramping_hist, digits=1))
        println()
        println("Maximum Production - Model (MW): ", round.(maxprod, digits=2))
        println("Maximum Production - Historical Data (MW): ", round(maximum(hist_prod_tot), digits=2))
        println("Minimum Production - Model (MW): ", round.(minprod, digits=2))
        println("Minimum Production - Historical Data (MW): ", round(minimum(hist_prod_tot), digits=2))
        println()
        println("Prod. above 85% of total capacity - Model (nr of h): ", highnr) 
        println("Prod. above 85% of total capacity - Historical Data (nr of h): ", count(x -> x >= 0.85*installed_capacity, hist_prod_tot)) 
        println("Prod. below 15% of total capacity - Model (nr of h): ", lownr) 
        println("Prod. below 15% of total capacity - Historical Data (nr of h): ", count(x -> x <= 0.15*installed_capacity, hist_prod_tot)) 
        println()
        println("RMSD on production, hourly values: ", round.(RMSD, digits=2))
        #println("RMSD on head, hourly values: ", round.(RMSD_head, digits=2))
    end

    letters = [layouts[name][1] for name in allrunnames]
    namedcolors = [layouts[name][2] for name in allrunnames]
    fillalphas = [layouts[name][5] for name in allrunnames]
    groupedruns = repeat(letters, outer=length(years)) |> CategoricalArray # same hack again
    levels!(groupedruns, letters)
    
    groupedbar(repeat(string.(years), inner=length(allrunnames)), RMSE; hover=RMSE,
            group = groupedruns, xlabel = "Year", ylabel = "RMSD", 
            #title = "RMSD on hourly production - model vs hist. data", 
            opt..., size=(450*2+50*2*length(years), 500*2),
            fillalpha = repeat(fillalphas, outer=length(years)),
            color = repeat(namedcolors, outer=length(years))) |> display
end
