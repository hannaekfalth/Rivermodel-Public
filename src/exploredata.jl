function exploredata(; plants=1:17, fit=0.7, span=0.25, lag=0, plottype=:prodeq, variant=:linear,
            years=[2008,2015], ylim=(0.5,1.0), plotsize=(1000,600), numsegments=2, showplots=true,
            adjusteta=false, nlags=3, ndownstreamlags=1, scaleeta=false, etaquantiles=nothing)
    vars = load("$DATAFOLDER/skelleftedata_raw.jld2")
    @unpack hours, production, forebay_level, tail_level, discharge, reservoir_weekly, inflow_daily, lake_levels = vars

    year = Dates.year.(hours)

    up = @view forebay_level[:,16]
    quick_clean!(up, mean(up[up.>0]))
    low = @view tail_level[:,16]
    quick_clean!(low, mean(low[low.>0]))
    
    head = forebay_level - tail_level
    calc_elec = 9.81*997*1e-6 * discharge .* head
    eta = production ./ calc_elec

    PLANT = riverdata(:Skellefte)[1]
    plantnames = [p.name for p in PLANT]
    SEGMENT = 1:numsegments
    nplants = length(PLANT)
    eta_seg, dc_seg =
        if etaquantiles === nothing
            estimate_eta_discharge(vars, plantnames, SEGMENT, PLANT; scaleeta)
        else
            estimate_eta_discharge(vars, plantnames, SEGMENT, PLANT, etaquantiles; scaleeta)
        end

    span0 = span
    lagdischarge = lags(discharge, 1)
    # lagdischarge2 = lags(discharge, 2)
    okyears = [y in years for y in year]

    meanabserror = zeros(nplants)

    for i in plants
        f = okyears .& (discharge[:, i] .>= 0)
        span = span0/(sum(f)/8760)
        if plottype == :prodeq
            maxelec = maximum(calc_elec[f, i])
            ff = f .& (calc_elec[:, i] .< fit * maxelec)
            coeff = linreg(calc_elec[ff, i], production[ff, i])  # or just [ones(length(c)) c] \ p
            plot_prodeq(calc_elec[f, i], production[f, i], PLANT[i], coeff; size)
            x = LinRange(0, maxelec, 2)
            y = coeff[1] .+ coeff[2]*x
            plot!(x, y, c=:black)
            plot_loessfit!(calc_elec[f, i], production[f, i]; span) |> display
        elseif plottype == :discharge_prod
            plot_eta(discharge[f, i], production[f, i], PLANT[i]; plotsize, ylim=:auto, loess=false, etafilter=false)
            md, cap = PLANT[i].maxdischarge, PLANT[i].capacity
            mxd, mxp = maximum(discharge[f, i]), maximum(production[f, i])
            plot!([md, md], [0, cap], c=2)
            plot!([0, md], [cap, cap], c=2)
            plot!([mxd, mxd], [0, mxp], c=3)
            plot!([0, mxd], [mxp, mxp], c=3)
            mxd, mxp = quantile(discharge[f, i], 0.999), quantile(production[f, i], 0.999)
            plot!([mxd, mxd], [0, mxp], c=3)
            plot!([0, mxd], [mxp, mxp], c=3) |> display
        elseif plottype == :discharge_proderror
            # assume forebay_level known, estimate tail_level, head & production
            # compare constant head, taylorhead & bilinear head
            mxd, mxp = maximum(discharge[f, i]), maximum(production[f, i])
            dd = 0:0.01:mxd
            coeff = [ones(sum(f)) discharge[f, i]] \ tail_level[f, i]
            predlow = coeff[1] .+ coeff[2]*dd
            coeff = [ones(sum(f)) discharge[f, i] discharge[f, i].^2 discharge[f, i].^3 discharge[f, i].^4] \ tail_level[f, i]
            predlow_poly = coeff[1] .+ coeff[2]*dd.+ coeff[3]*dd.^2 .+ coeff[4]*dd.^3 .+ coeff[5]*dd.^4
            head = forebay_level - tail_level
            calcprod = 9.81*997*1e-6 * discharge .* head
            eta = production ./ calc_elec
            plot_eta(discharge[f, i], production[f, i], PLANT[i]; plotsize, ylim=:auto, loess=false, etafilter=false)
            md, cap = PLANT[i].maxdischarge, PLANT[i].capacity
            
            plot!([md, md], [0, cap], c=2)
            plot!([0, md], [cap, cap], c=2)
            plot!([mxd, mxd], [0, mxp], c=3)
            plot!([0, mxd], [mxp, mxp], c=3)
            mxd, mxp = quantile(discharge[f, i], 0.999), quantile(production[f, i], 0.999)
            plot!([mxd, mxd], [0, mxp], c=3)
            p = plot!([0, mxd], [mxp, mxp], c=3)
            segmax = cumsum(dc_seg[i,:])
            segmax1 = [0; segmax[1:end-1]]

            dds = [d .> segmax[s] ? segmax[s]-segmax1[s] : max(0, d-segmax1[s]) for d in dd, s = 1:numsegments]
            pp = sumdrop(9.81*997*1e-6 * dds .* (mean(forebay_level[f, i]) .- predlow) .* eta_seg[i,:]', dims=2)
            ppp = sumdrop(9.81*997*1e-6 * dds .* (mean(forebay_level[f, i]) .- predlow_poly) .* eta_seg[i,:]', dims=2)
            pph = sumdrop(9.81*997*1e-6 * dds .* mean(head[f, i]) .* eta_seg[i,:]', dims=2)
            plot!(dd, [pph pp ppp], c=[:red :black :purple], ylim=ylims(p)) |> display
        elseif plottype == :head_eta
            f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0)
            low, high = quantile(head[f, i], [0.001, 0.999])
            # dc_high = quantile(discharge[f, i], 0.9)
            f = f .& (low .< head[:, i] .< high) # .& (discharge[:, i] .> dc_high)
            plot_eta(head[f, i], eta[f, i], PLANT[i]; ylim, plotsize, span, loess=false) |> display
        elseif plottype == :head_prod
            f = okyears .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0)
            low, high = quantile(head[f, i], [0.001, 0.999])
            f = f .& (low .< head[:, i] .< high)
            plot_eta(head[f, i], production[f, i], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false) |> display
        elseif plottype == :discharge_eta
            plot_eta(discharge[f, i], eta[f, i], PLANT[i]; ylim, plotsize, span)
            cumdc = cumsum(dc_seg[i,:])
            weightedeta = cumsum(dc_seg[i,:] .* eta_seg[i,:]) ./ cumdc
            xx = repeat([0; cumdc], inner=2)[2:end-1]
            yy = repeat(weightedeta, inner=2)
            plot!(xx, yy, c=:black) |> display
        elseif plottype == :discharge_effdischarge
            plot_eta(discharge[f, i], discharge[f, i].*eta[f, i], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false)
            mxd = maximum(discharge[f, i])
            dd = 0:0.01:mxd
            segmax = cumsum(dc_seg[i,:])
            segmax1 = [0; segmax[1:end-1]]          
            dds = [d .> segmax[s] ? segmax[s]-segmax1[s] : max(0, d-segmax1[s]) for d in dd, s = 1:numsegments]
            effdds = dds .* eta_seg[i,:]'
            plot!(dd, [dd*eta_seg[i,1] sumdrop(effdds, dims=2)], c=[:red :black]) |> display
        elseif plottype == :time_head
            f = okyears .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0)
            low, high = quantile(head[f, i], [0.001, 0.999])
            plot_eta(hours[f], head[f, i], PLANT[i]; ylim=(low*0.95, high*1.05), plotsize, span, loess=false, etafilter=false,
                        titlesuffix="$(round(high-low, digits=1)) m")
            plot!([hours[f][1], hours[f][end]], [low, low], c=2)
            plot!([hours[f][1], hours[f][end]], [high, high], c=2) |> display
        elseif plottype == :time_prod
            f = okyears
            fff = okyears .& ((year .== 2008) .| (year .== 2015))
            cap = PLANT[i].capacity
            mxp = maximum(production[fff, i])
            mx_pseg = 9.81*997*1e-6 * dc_seg[i, 1] .* PLANT[i].reportedhead * eta_seg[i, 1]
            plot_eta(hours[f], production[f, i], PLANT[i]; ylim=(0, max(mxp,cap)*1.05), plotsize, span, loess=false, etafilter=false)
            plot!([hours[f][1], hours[f][end]], [cap, cap], c=2)
            plot!([hours[f][1], hours[f][end]], [mxp, mxp], c=3)
            plot!([hours[f][1], hours[f][end]], [mx_pseg, mx_pseg], c=:black) |> display
        elseif plottype == :time_calcdischarge
            f = okyears
            fff = okyears .& ((year .== 2008) .| (year .== 2015))
            calcdischarge = production[f, i] ./ head[f, i] / (9.81*997*1e-6) / eta_seg[i, 1]
            md = PLANT[i].maxdischarge
            mxd = maximum(discharge[fff, i])
            mx_dcseg = dc_seg[i, 1]
            p = plot_eta(hours[f], calcdischarge, PLANT[i]; ylim=(0, max(mxd,md)*1.05), plotsize, span, loess=false, etafilter=false)
            plot!([hours[f][1], hours[f][end]], [md, md], c=2)
            plot!([hours[f][1], hours[f][end]], [mxd, mxd], c=3) 
            plot!([hours[f][1], hours[f][end]], [mx_dcseg, mx_dcseg], c=:black) |> display
        elseif plottype == :time_calcdischarge_alt
            f = okyears .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0)
            calcdischarge = production[f, i] ./ head[f, i] / (9.81*997*1e-6) / eta_seg[i, 1]
            maxhead = quantile(head[f, i], 0.999)
            etafactor = (1 .- 0.24*maxhead./head[f, i]) / (1 - 0.24)
            # meanhead = mean(head[f, i])
            # etafactor = (1 .- 0.27*meanhead./head[f, i]) / (1 - 0.27)
            mx = quantile(calcdischarge./etafactor, 0.999)
            calcdischarge = adjusteta ? calcdischarge./etafactor : calcdischarge
            plot_eta(hours[f], calcdischarge, PLANT[i]; ylim=(0, mx*1.05), plotsize, span, loess=false, etafilter=false) |> display
        elseif plottype == :time_discharge
            f = okyears
            fff = okyears .& ((year .== 2008) .| (year .== 2015))
            md = PLANT[i].maxdischarge
            mxd = maximum(discharge[fff, i])
            mx_dcseg = dc_seg[i, 1]
            p = plot_eta(hours[f], discharge[f, i], PLANT[i]; ylim=(0, max(mxd,md)*1.05), plotsize, span, loess=false, etafilter=false)
            plot!([hours[f][1], hours[f][end]], [md, md], c=2)
            plot!([hours[f][1], hours[f][end]], [mxd, mxd], c=3) 
            plot!([hours[f][1], hours[f][end]], [mx_dcseg, mx_dcseg], c=:black) |> display
        elseif plottype == :time_lower
            f = okyears
            plot_eta(hours[f], tail_level[f, i], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false)
            f = okyears .& (discharge[:, i] .> 0)
            if variant == :linear
                coeff = [ones(sum(f)) discharge[f, i]] \ tail_level[f, i]
                pred = coeff[1] .+ coeff[2]*discharge[f, i]
                scatter!(hours[f], pred, markersize=1, markerstrokewidth=0) |> display
            elseif variant == :linearupper
                coeff = [ones(sum(f)) discharge[f, i] forebay_level[f, i]] \ tail_level[f, i]
                pred = coeff[1] .+ coeff[2]*discharge[f, i] .+ coeff[3]*forebay_level[f, i]
                scatter!(hours[f], pred, markersize=1, markerstrokewidth=0) |> display
            elseif variant == :quadratic
                coeff = [ones(sum(f)) discharge[f, i] discharge[f, i].^2] \ tail_level[f, i]
                pred = coeff[1] .+ coeff[2]*discharge[f, i] .+ coeff[3]*discharge[f, i].^2
                scatter!(hours[f], pred, markersize=1, markerstrokewidth=0) |> display
            elseif variant == :poly4
                coeff = [ones(sum(f)) discharge[f, i] discharge[f, i].^2 discharge[f, i].^3 discharge[f, i].^4] \ tail_level[f, i]
                pred = coeff[1] .+ coeff[2]*discharge[f, i] .+ coeff[3]*discharge[f, i].^2 .+ coeff[4]*discharge[f, i].^3 .+ coeff[5]*discharge[f, i].^4
                scatter!(hours[f], pred, markersize=1, markerstrokewidth=0) |> display
            end
        elseif plottype == :time_upper
            f = okyears
            plot_eta(hours[f], forebay_level[f, i], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false)
            limits = [PLANT[i].forebayhigh PLANT[i].forebaylow]
            plot!([extrema(hours[f])...], [limits; limits], c=:green) |> display
        elseif plottype == :lagdischarge_discharge
            f = f .& (lagdischarge[:, i] .> 0)
            plot_eta(lagdischarge[f, i], discharge[f, i], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false) |> display
        elseif plottype == :dischargelagdiff_lower
            f = f .& (tail_level[:, i] .> 0) .& (lagdischarge[:, i] .> 0)
            dischargelagdiff = discharge[:,i] - lagdischarge[:,i]
            plot_eta(dischargelagdiff[f], tail_level[f, i], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false) |> display
        elseif plottype == :discharge_lower
            i == 17 && continue
            f = f .& (tail_level[:, i] .> 0)
            maxhead = quantile(head[f, i], 0.999)
            etafactor = (1 .- 0.24*maxhead./head[:, i]) / (1 - 0.24)
            dc = adjusteta && i != 3 ? discharge[:, i] ./ etafactor : discharge[:, i]
            lagdc = hcat([lags(dc, lag) for lag = 1:nlags]...)
            plot_eta(dc[f], tail_level[f, i], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false)
            for iter = 1:2
                if iter == 2
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    f = f .& (forebay_level[:,downstream] .> 0) .& (tail_level[:, i] .>= forebay_level[:,downstream])
                end
                if variant == :linear
                    coeff = [ones(sum(f)) dc[f]] \ tail_level[f, i]
                    downstream = i == 1 ? 3 : i+1
                    # display([tail_level[:, 3] forebay_level[:,3]])
                    x = LinRange(extrema(dc[f])..., 500)
                    y = coeff[1] .+ coeff[2]*x
                    plot!(x, y)
                    pred = coeff[1] .+ coeff[2]*dc[f]
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper
                    f = f .& (forebay_level[:, i] .> 0)
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] .+ coeff[3]*forebay_level[f, i]
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_downstreamupper
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0)
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, downstream]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] .+ coeff[3]*forebay_level[f, downstream]
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_ar1discharge
                    ar1_dc = ar1(dc[f], 0.6)
                    coeff = [ones(sum(f)) ar1_dc] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*ar1_dc
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_delta_downstreamupper
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    downstreamupper = forebay_level[:, downstream]
                    downstreamupper2 = [0; downstreamupper[1:end-1]]
                    resdiff = [0; diff(downstreamupper)]
                    f = f .& (downstreamupper .> 0) .& (downstreamupper2 .> 0) .& isfinite.(resdiff)
                    ar1_resdiff = ar1(resdiff[f], 0.6)
                    coeff = [ones(sum(f)) dc[f] ar1_resdiff] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] .+ coeff[3]*ar1_resdiff
                    # scatter!(dcok, pred, markersize=1, markerstrokewidth=0)
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper_downstreamupper
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0)
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i] forebay_level[f, downstream]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*forebay_level[f, i] + coeff[4]*forebay_level[f, downstream]
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper_downstreamupper_quad
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0)
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i] forebay_level[f, downstream] dc[f].^2] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] .+ coeff[3]*forebay_level[f, i] .+ coeff[4]*forebay_level[f, downstream] .+ coeff[5]*dc[f].^2
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper_downstreamupper_lag
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0) .&
                                vec(prod(lagdc .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i] forebay_level[f, downstream] lagdc[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*forebay_level[f, i] + coeff[4]*forebay_level[f, downstream] + lagdc[f, :] * coeff[5:4+nlags]
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_downstreamupper_lag
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0) .&
                                vec(prod(lagdc .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, downstream] lagdc[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*forebay_level[f, downstream] + lagdc[f, :] * coeff[4:3+nlags]
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper_downstreamupper_lag_dlag
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    lagdownstream = hcat([lags(forebay_level[:, downstream], lag) for lag = 1:ndownstreamlags]...)
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0) .&
                                vec(prod(lagdc .> 0, dims=2)) .& vec(prod(lagdownstream .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i] forebay_level[f, downstream] lagdc[f, :] lagdownstream[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*forebay_level[f, i] + coeff[4]*forebay_level[f, downstream] +
                                lagdc[f, :] * coeff[5:4+nlags] + lagdownstream[f, :] * coeff[5+nlags:4+nlags+ndownstreamlags] 
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper_downstreamupper_quadlag
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0) .&
                                vec(prod(lagdc .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i] forebay_level[f, downstream] dc[f].^2 lagdc[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*forebay_level[f, i] + coeff[4]*forebay_level[f, downstream] +
                                coeff[5]*dc[f].^2 + lagdc[f, :] * coeff[6:5+nlags] 
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper_downstreamupper_quadlag_dlag
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    lagdownstream = hcat([lags(forebay_level[:, downstream], lag) for lag = 1:ndownstreamlags]...)
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0) .&
                                vec(prod(lagdc .> 0, dims=2)) .& vec(prod(lagdownstream .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i] forebay_level[f, downstream] dc[f].^2 lagdc[f, :] lagdownstream[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*forebay_level[f, i] + coeff[4]*forebay_level[f, downstream] +
                                coeff[5]*dc[f].^2 + lagdc[f, :] * coeff[6:5+nlags] + lagdownstream[f, :] * coeff[6+nlags:5+nlags+ndownstreamlags] 
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linear_upper_downstreamupper_dlag
                    i == 17 && continue
                    downstream = i == 1 ? 3 : i+1
                    lagdownstream = hcat([lags(forebay_level[:, downstream], lag) for lag = 1:ndownstreamlags]...)
                    f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0) .&
                                vec(prod(lagdownstream .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] forebay_level[f, i] forebay_level[f, downstream] lagdownstream[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*forebay_level[f, i] + coeff[4]*forebay_level[f, downstream] +
                                lagdownstream[f, :] * coeff[5:4+ndownstreamlags] 
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    # @show coeff
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :quadratic
                    coeff = [ones(sum(f)) dc[f] dc[f].^2] \ tail_level[f, i]
                    x = LinRange(extrema(dc[f])..., 500)
                    y = coeff[1] .+ coeff[2]*x .+ coeff[3]*x.^2
                    plot!(x, y)
                    pred = coeff[1] .+ coeff[2]*dc[f] .+ coeff[3]*dc[f].^2
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :linearquadratic
                    coeff1 = [ones(sum(f)) dc[f]] \ tail_level[f, i]
                    x = LinRange(extrema(dc[f])..., 500)
                    y1 = coeff1[1] .+ coeff1[2]*x
                    coeff2 = [ones(sum(f)) dc[f] dc[f].^2] \ tail_level[f, i]
                    y2 = coeff2[1] .+ coeff2[2]*x .+ coeff2[3]*x.^2
                    plot!(x, [y1 y2])
                elseif variant == :linearlag
                    f = f .& vec(prod(lagdc .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] lagdc[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + lagdc[f, :] * coeff[3:end]
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :quadraticlag
                    f = f .& (lagdischarge[:, i] .> 0)
                    coeff = [ones(sum(f)) dc[f] lagdc[f] dc[f].^2 lagdc[f].^2] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] .+ coeff[3]*lagdc[f] .+ coeff[4]*dc[f].^2 .+ coeff[5]*lagdc[f].^2
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :poly4
                    coeff = [ones(sum(f)) dc[f] dc[f].^2 dc[f].^3 dc[f].^4] \ tail_level[f, i]
                    x = LinRange(extrema(dc[f])..., 500)
                    y = coeff[1] .+ coeff[2]*x .+ coeff[3]*x.^2 .+ coeff[4]*x.^3 .+ coeff[5]*x.^4
                    plot!(x, y)
                    pred = coeff[1] .+ coeff[2]*dc[f] .+ coeff[3]*dc[f].^2 .+ coeff[4]*dc[f].^3 .+ coeff[5]*dc[f].^4
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                elseif variant == :poly4lag
                    f = f .& vec(prod(lagdc .> 0, dims=2))
                    coeff = [ones(sum(f)) dc[f] dc[f].^2 dc[f].^3 dc[f].^4 lagdc[f, :]] \ tail_level[f, i]
                    pred = coeff[1] .+ coeff[2]*dc[f] + coeff[3]*dc[f].^2 + coeff[4]*dc[f].^3 + coeff[5]*dc[f].^4 + lagdc[f, :] * coeff[6:end]
                    scatter!(dc[f], pred, markersize=1, markerstrokewidth=0)
                    meanabserror[i] = mean(abs.(tail_level[f, i] .- pred))
                end
            end
            showplots && display(plot!())
        elseif plottype == :discharge_lower_head
            f = f .& (tail_level[:, i] .> 0)
            colorlimits = Tuple(quantile(head[f, i], [0.1, 0.9]))
            plot_color2d(discharge[f, i], tail_level[f, i], head[f, i], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :discharge_lower_upper
            f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0)
            colorlimits = Tuple(quantile(forebay_level[f, i], [0.1, 0.9]))
            plot_color2d(discharge[f, i], tail_level[f, i], forebay_level[f, i], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :discharge_lower_downstreamupper
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0)
            colorlimits = Tuple(quantile(forebay_level[f,downstream], [0.1, 0.9]))
            plot_color2d(discharge[f, i], tail_level[f, i], forebay_level[f,downstream], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :downstreamupper_lower_discharge
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:,downstream] .> 0)
            colorlimits = Tuple(quantile(discharge[f, i], [0.1, 0.9]))
            plot_color2d(forebay_level[f,downstream], tail_level[f, i], discharge[f, i], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :discharge_lower_dischargeupstream
            i <= 3 && continue
            dischargeupstream = discharge[:,i-1]
            f = f .& (tail_level[:, i] .> 0) .& (dischargeupstream .> 0)
            colorlimits = Tuple(quantile(dischargeupstream[f], [0.1, 0.9]))
            plot_color2d(discharge[f, i], tail_level[f, i], dischargeupstream[f], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :upper_lower_discharge
            f = f .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0)
            colorlimits = Tuple(quantile(discharge[f, i], [0.1, 0.9]))
            plot_color2d(forebay_level[f, i], tail_level[f, i], discharge[f, i], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :discharge_lower_eta
            f = f .& (tail_level[:, i] .> 0)
            colorlimits = Tuple(quantile(eta[f, i], [0.1, 0.9]))
            plot_color2d(discharge[f, i], tail_level[f, i], eta[f, i], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :discharge_lower_day
            f = f .& (tail_level[:, i] .> 0)
            day = dayofyear.(hours)
            colorlimits = Tuple(quantile(day[f], [0.1, 0.9]))
            plot_color2d(discharge[f, i], tail_level[f, i], day[f], PLANT[i]; plotsize, clims=colorlimits) |> display
        elseif plottype == :autocor_discharge
            ac = autocor(discharge[f, i], collect(0:24))
            bar(ac) |> display
        elseif plottype == :pacf_discharge
            pac = pacf(discharge[f, i], collect(0:24))
            bar(pac) |> display
        elseif plottype == :crosscor_discharge_lower
            cc = crosscor(discharge[f, i], tail_level[f, i], collect(-24:24))
            bar(cc) |> display
        elseif plottype == :autocor_diffdischarge
            ac = autocor(diff(discharge[f, i]), collect(0:24))
            bar(ac) |> display
        elseif plottype == :pacf_diffdischarge
            pac = pacf(diff(discharge[f, i]), collect(0:24))
            bar(pac) |> display
        elseif plottype == :prod_downstreamupper
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            f = okyears .& (production[:, i] .> 0) .& (forebay_level[:, downstream] .> 0)
            plot_eta(production[f, i], forebay_level[f, downstream], PLANT[i]; ylim=:auto, plotsize, span, loess=false, etafilter=false) |> display
        elseif plottype == :crosscor_prod_downstreamupper
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            vprod = copy(production[okyears, i])
            downstreamupper = forebay_level[okyears, downstream]
            f1 = (vprod .>= 0) .& isfinite.(vprod)
            f2 = (downstreamupper .> 0) .& isfinite.(downstreamupper)
            # println(sum(f),"   ",sum(okyears))
            replacevalues!(vprod, f1)  # replace bad values with mean data so lags don't get messed up
            replacevalues!(downstreamupper, f2)
            lag = collect(-72:72)
            cc = crosscor(vprod, downstreamupper, lag)
            bar(lag, cc, title="$(PLANT[i].name)-$(PLANT[downstream].name)", xticks=-70:10:70) |> display
        elseif plottype == :crosscor_prod_delta_downstreamupper
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            vprod = copy(production[okyears, i])
            downstreamupper = forebay_level[okyears, downstream]
            downstreamupper2 = [0; downstreamupper[1:end-1]]
            resdiff = [0; diff(downstreamupper)]
            f1 = (vprod .>= 0) .& isfinite.(vprod)
            f2 = (downstreamupper .> 0) .& (downstreamupper2 .> 0) .& isfinite.(resdiff)
            # println(sum(f),"   ",sum(okyears))
            replacevalues!(vprod, f1)  # replace bad values with mean data so lags don't get messed up
            replacevalues!(resdiff, f2)
            lag = collect(-72:72)
            cc = crosscor(vprod, resdiff, lag)
            corr = max.(0, round.(Int, cc[73:73+30] * 10000))  # lags 0-30
            @show corr
            bar(lag, cc, title="$(PLANT[i].name)-$(PLANT[downstream].name)", xticks=-70:10:70) |> display
        elseif plottype == :crosscor_prod_delta_lower
            vprod = copy(production[okyears, i])
            lower = tail_level[okyears, i]
            lower2 = [0; lower[1:end-1]]
            resdiff = [0; diff(lower)]
            f1 = (vprod .>= 0) .& isfinite.(vprod)
            f2 = (lower .> 0) .& (lower2 .> 0) .& isfinite.(resdiff)
            # println(sum(f),"   ",sum(okyears))
            replacevalues!(vprod, f1)  # replace bad values with mean data so lags don't get messed up
            replacevalues!(resdiff, f2)
            lag = collect(-72:72)
            cc = crosscor(vprod, resdiff, lag)
            bar(lag, cc, title="$(PLANT[i].name)", xticks=-70:10:70) |> display
        elseif plottype == :crosscor_prod_lower
            vprod = copy(production[okyears, i])
            lower = tail_level[okyears, i]
            f1 = (vprod .>= 0) .& isfinite.(vprod)
            f2 = (lower .> 0) .& isfinite.(lower)
            # println(sum(f),"   ",sum(okyears))
            replacevalues!(vprod, f1)  # replace bad values with mean data so lags don't get messed up
            replacevalues!(lower, f2)
            lag = collect(-72:72)
            cc = crosscor(vprod, lower, lag)
            bar(lag, cc, title="$(PLANT[i].name)", xticks=-70:10:70) |> display
        elseif plottype == :discharge_downstreameta
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            downstreamhead = forebay_level[f, i] - forebay_level[f, downstream]
            downstreamelec = 9.81*997*1e-6 * discharge[f, i] .* downstreamhead
            downstreameta = production[f, i] ./ downstreamelec
            plot_eta(discharge[f, i], downstreameta, PLANT[i]; ylim, plotsize, span) |> display
        elseif plottype == :discharge_head_eta
            plot_color2d(discharge[f, i], head[f, i], eta[f, i], PLANT[i]; size) |> display
        elseif plottype == :discharge_head
            plot_eta(discharge[f, i], head[f, i], PLANT[i]; plotsize, ylim=:auto, loess=false, etafilter=false)
            plot!([extrema(discharge[f, i])...], [1,1] * [PLANT[i].reportedhead mean(head[f, i])]) |> display
        elseif plottype == :discharge_downstreamhead
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            downstreamhead = forebay_level[f, i] - forebay_level[f, downstream]
            plot_eta(discharge[f, i], downstreamhead, PLANT[i]; plotsize, ylim=:auto, loess=false, etafilter=false) |> display
        elseif plottype == :discharge_downstreamdiff
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            downstreamdiff = tail_level[f, i] - forebay_level[f, downstream]
            plot_eta(discharge[f, i], downstreamdiff, PLANT[i]; plotsize, ylim=:auto, loess=false, etafilter=false) |> display
        elseif plottype == :discharge_head_downstreamdiff
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            downstreamdiff = tail_level[f, i] - forebay_level[f, downstream]
            colorlimits = Tuple(quantile(downstreamdiff, [0.1, 0.9]))
            try
                plot_color2d(discharge[f, i], head[f, i], downstreamdiff, PLANT[i]; plotsize, clims=colorlimits) |> display
            catch
                @warn "Error on plot $i"
            end
        elseif plottype == :lower_downstreamdiff
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            f = okyears .& (tail_level[:, i] .> 0) .& (forebay_level[:, downstream] .> 0)
            downstreamdiff = tail_level[f, i] - forebay_level[f, downstream]
            plot_eta(tail_level[f, i].-minimum(tail_level[f, i]), downstreamdiff, PLANT[i];
                        ylim=:auto, plotsize, span, loess=false, etafilter=false) |> display
        elseif plottype == :time_downstreamdiff
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            f = okyears .& (tail_level[:, i] .> 0) .& (forebay_level[:, downstream] .> 0)
            downstreamdiff = tail_level[f, i] - forebay_level[f, downstream]
            plot_eta(hours[f], downstreamdiff, PLANT[i];
                        ylim=:auto, plotsize, span, loess=false, etafilter=false) |> display
        elseif plottype == :lower_downstreamupper
            i == 17 && continue
            downstream = i == 1 ? 3 : i+1
            f = okyears .& (tail_level[:, i] .> 0) .& (forebay_level[:, downstream] .> 0)
            p = plot_eta(tail_level[f, i], forebay_level[f, downstream], PLANT[i];
                        ylim=:auto, plotsize, span, loess=false, etafilter=false)
            xlimits = quantile(tail_level[f, i], [0.001, 0.999])
            ylimits = quantile(forebay_level[f, downstream], [0.001, 0.999])
            plot!(xlimits, xlimits; xlim=xlimits, ylim=ylimits) |> display
        elseif plottype == :time_levels
            f = okyears .& (tail_level[:, i] .> 0) .& (forebay_level[:, i] .> 0)
            plant = PLANT[i]
            title = string(plant.name) * " $(plant.capacity) MW, $(plant.maxdischarge) m3/s"
            scatter(hours[f], [forebay_level[f, i] tail_level[f, i]];
                        markersize=1, markerstrokewidth=0, legend=:none, alpha=1, 
                        tickfont=12*plotsize[1]÷1000, ylim=:auto, plotsize, title)
            if i <= 5
                plot!(hours[okyears], lake_levels[okyears, i];
                            markersize=1, markerstrokewidth=0, legend=:none, alpha=1, c=:black)
            end
            !any(f) && (display(plot!()); continue)
            lower = quantile(tail_level[f, i], [0.001, 0.999])
            upper = quantile(forebay_level[f, i], [0.001, 0.999])
            plot!([extrema(hours[f])...], [lower[1], lower[1]], c=:green)
            plot!([extrema(hours[f])...], [lower[2], lower[2]], c=:green)
            plot!([extrema(hours[f])...], [upper[1], upper[1]], c=:green)
            plot!([extrema(hours[f])...], [upper[2], upper[2]], c=:green, ticks=:native) |> display
        else
            @warn "Unrecognized plot type."
        end
    end
    if plottype == :discharge_lower
        meanabserror = round.(meanabserror, digits=3)
        any(meanabserror .!= 0) && @show meanabserror
        @show mean(meanabserror[[1; 2; 4:end-1]])
    end
    nothing
end

function plot_color2d(discharge, head, eta, plant; plotsize=(750,400), clims=(0.7,0.95))
    title = string(plant.name) * " $(plant.capacity) MW, $(plant.maxdischarge) m3/s"
    opt = (markersize=1, markerstrokewidth=0, legend=:none, alpha=1.0, colorbar=true, title, size)
    scatter(discharge, head; opt..., marker_z=eta, color=palette(:rainbow,5), clims=clims)
end

function plot_eta(x, eta, plant; ylim=(0.5,1.0), plotsize=(750,400), span=0.25, loess=true, etafilter=true, titlesuffix="")
    title = string(plant.name) * " $(plant.capacity) MW, $(plant.maxdischarge) m3/s $titlesuffix"
    opt = (markersize=1, markerstrokewidth=0, legend=:none, ticks=:native, alpha=0.3, tickfont=12*plotsize[1]÷1000, title, size=plotsize, ylim)
    f = etafilter ? eta .< 1 : fill(true, length(eta))
    p = scatter(x[f], eta[f]; opt...)
    loess && plot_loessfit!(x[f], eta[f]; span=span)
    return p
end

function plot_prodeq(calc_elec, production, plant, coeff; plotsize=(750,400))
    stats = stringify(η=coeff[2], Δx=-coeff[1]/coeff[2], R²=r_squared(production,calc_elec))
    title = string(plant.name) * " $(plant.capacity) MW, $(plant.maxdischarge) m3/s ($stats)"
    opt = (markersize=1, markerstrokewidth=0, legend=:none, alpha=0.4, title, size)
    maxelec = maximum(calc_elec)
    scatter(calc_elec, production; opt...)
    # plot!([0,maxelec], [0,maxelec], c=2)
    # plot!([0,maxelec], [0,0.9*maxelec], c=2)
    h = fit(Histogram, calc_elec, nbins=200)
    hist = Histogram(h.edges,maxelec/2*h.weights/maximum(h.weights),h.closed,h.isdensity)
    plot!(hist, c=:black, alpha=0.15)
end
