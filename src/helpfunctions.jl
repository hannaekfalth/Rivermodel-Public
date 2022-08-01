using LinearAlgebra, Loess, StatsPlots

shift(v, t) = v[mod1(t - v[1] + 1, length(v))]
sumdrop(x; dims) = (dims > ndims(x)) ? x : dropdims(sum(x; dims); dims)
meandrop(x; dims) = (dims > ndims(x)) ? x : dropdims(mean(x; dims); dims)
lookup(dict, keys::AbstractArray) = getindex.(Ref(dict), keys)

hourofyear(dt) = (dayofyear(dt) - 1)*24 + 1 + hour(dt)
timerange(year, start, stop) =
    hourofyear(DateTime("$(year)-$(start)T00")):hourofyear(DateTime("$(year)-$(stop)T23"))

function moving_average(array, n)
    output=zeros(length(array))

    for i=1:length(array)
        output[i] = sum(shift(array, i) for i in i:i+n-1)/n
    end

    return output
end

# https://discourse.julialang.org/t/efficient-way-of-doing-linear-regression/31232/28
function linreg(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:AbstractFloat}
    (N = length(x)) == length(y) || throw(DimensionMismatch())
    ldiv!(cholesky!(Symmetric([T(N) sum(x); zero(T) sum(abs2, x)], :U)), [sum(y), dot(x, y)])
end

r_squared(y, yfit) = 1 - sum((y .- yfit).^2) / sum((y .- mean(y)).^2)

function plot_loessfit!(x, y; span=0.25)
    model = loess(x, y, span=span)
    xx = LinRange(extrema(x)..., 1000)
    yfit = Loess.predict(model, xx)
    plot!(xx, yfit, c=:red)
end

# for pretty printing a NamedTuple
stringify(; sigdigits=3, args...) = string(map(x -> round(x, sigdigits=sigdigits), NamedTuple(args)))

normalize01(x) = (ex = extrema(x); return (x .- ex[1]) / (ex[2] - ex[1]))
lags(x,n) = [zeros(eltype(x), n, size(x)[2:end]...); selectdim(x, 1, 1:size(x,1)-n)]
neglags(x,n) = [x[1+n:end, :]; zeros(n, size(x,2))]

replacevalues!(x, f) = (x[.!f] .= mean(x[f]))

function ar1(x, rho)
    y = similar(x)
    y[1] = x[1]
    for i=2:length(x)
        y[i] = rho*y[i-1] + (1-rho)*x[i]
    end
    return y / sqrt((1-rho)^2 / (1-rho^2))
end

# positivequantile(x::AbstractVector, q::Real) = positivequantile(x, q)[1] 
function positivequantile(x::AbstractVector, q)
    f = (x .> 0)
    return any(f) ? quantile(x[f], q) : (q isa Vector) ? fill(-Inf, length(q)) : -Inf
end

function quick_clean!(vect, avg; forbid_zero=true)
    last_ok = avg
    for (i, e) in enumerate(vect)
        if e < 0 || (forbid_zero && e == 0)
            vect[i] = last_ok
        else
            last_ok = e
        end
    end
end

# test stuff...
# R.groupedbar(repeat(["2015","2019"], inner=3), [3 2 1; 6 5 4], label=["c" "b" "a"])
# ctg = R.CategoricalArray(["a", "b", "c", "a", "b", "c"])
# R.levels!(ctg, ["c", "b", "a"])
# R.groupedbar(repeat(["2015","2019"], inner=3), 1:6, group=ctg)
# R.groupedbar(repeat(["2015","2019"], inner=3), [1 2 3; 4 5 6]', group=ctg)

# errorbar demos
# years, names, data, low, high = R.testdata();
# R.bartest(years, names, data, low, high)
# R.bartest_CategoricalArray_v1(years, names, data, low, high)
# R.bartest_CategoricalArray_v2(years, names, data, low, high)

# years, names, data, low, high = R.testdata_flipped();
# R.bartest_flipped(years, names, data, low, high)
# R.bartest_CategoricalArray_v1_flipped(years, names, data, low, high)
# R.bartest_CategoricalArray_v2_flipped(years, names, data, low, high)

function testdata()
    years = 2015:4:2019
    names = ["C", "B", "A", "D", "E"]
    data = @aa [1000,2000] .+ 200*[1,2].*rand(length(years),length(names)), years, names
    display(data)
    low = (0.95 .- 0.1*rand(size(data)...)) .* data
    # low[:, 2:3] .= NaN
    high = (1.05 .+ 0.1*rand(size(data)...)) .* data
    return years, names, data, low, high
end

# bar chart with errorbars, simple version without CategoricalArray
function bartest(years, names, data, low, high; bar_width=0.8)
    plotly()
    y,n = size(data)
    groupedbar(repeat(string.(years), inner=n), data; bar_width, legend=:outertopright, label=permutedims(names))
    errorbars!(data, low, high; bar_width) |> display
end

# bar chart with errorbars using CategoricalArray
# works with standard errorbars!() but requires transposed data in groupedbar() call (but not in call to errorbars!)
function bartest_CategoricalArray_v1(years, names, data, low, high; bar_width=0.8)
    plotly()
    y,n = size(data)
    groupedruns = repeat(names, outer=y) |> CategoricalArray # same hack again
    levels!(groupedruns, names)
    groupedbar(repeat(string.(years), inner=n), data'; bar_width, legend=:outertopright, group=groupedruns)
    errorbars!(data, low, high; bar_width) |> display
end

# bar chart with errorbars using CategoricalArray when all data is naturally transposed
# requires an alternative version of errorbars!() adapted to transposed data
function bartest_CategoricalArray_v2(years, names, data, low, high; bar_width=0.8)
    plotly()
    y,n = size(data)
    groupedruns = repeat(names, outer=y) |> CategoricalArray # same hack again
    levels!(groupedruns, names)
    groupedbar(repeat(string.(years), inner=n), data'; bar_width, legend=:outertopright, group=groupedruns)
    errorbars_flipped!(data', low', high'; bar_width) |> display
end

function testdata_flipped()
    years = 2015:4:2019
    names = ["A", "B", "C", "D", "E"]
    data = @aa [1000 2000] .+ 200*[1 2].*rand(length(names),length(years)), names, years
    display(data)
    low = (0.95 .- 0.1*rand(size(data)...)) .* data
    # low[2:3, :] .= NaN
    high = (1.05 .+ 0.1*rand(size(data)...)) .* data
    return years, names, data, low, high
end

function bartest_flipped(years, names, data, low, high; bar_width=0.8)
    plotly()
    n,y = size(data)
    groupedbar(repeat(string.(years), inner=n), data'; bar_width, legend=:outertopright, label=permutedims(names))
    errorbars!(data', low', high'; bar_width) |> display
end

function bartest_CategoricalArray_v1_flipped(years, names, data, low, high; bar_width=0.8)
    plotly()
    n,y = size(data)
    groupedruns = repeat(names, outer=y) |> CategoricalArray # same hack again
    levels!(groupedruns, names)
    groupedbar(repeat(string.(years), inner=n), data; bar_width, legend=:outertopright, group=groupedruns)
    errorbars!(data', low', high'; bar_width) |> display
end

function bartest_CategoricalArray_v2_flipped(years, names, data, low, high; bar_width=0.8)
    plotly()
    n,y = size(data)
    groupedruns = repeat(names, outer=y) |> CategoricalArray # same hack again
    levels!(groupedruns, names)
    groupedbar(repeat(string.(years), inner=n), data; bar_width, legend=:outertopright, group=groupedruns)
    errorbars_flipped!(data, low, high; bar_width) |> display
end

# bar_width is total bar width of each group, so bar_width=1 will pack all bar groups together
function errorbars!(data, low, high; bar_width=0.8, line=(color=:black, linewidth=1),
                    label_low="warm start", marker_low=(:circle, :black, 2),
                    label_high="upper bound", marker_high=(:hline, :black, 4))
    r, c = size(data)
    xgroup = repeat(0:r-1, inner=c)
    barindex = repeat(0:c-1, outer=r)
    width = bar_width / c
    xfirstbar = (1 - bar_width)/2
    x = xgroup .+ (xfirstbar + width/2 .+ width*barindex)
    ylow = low'[:]
    yhigh = high'[:]
    xlines = [x x fill(NaN, r*c)]'[:]
    ylines = [data'[:] yhigh fill(NaN, r*c)]'[:]
    plot!(xlines, ylines; label="", line...)
    scatter!(x, yhigh, label=label_high, marker=marker_high)
    scatter!(x, ylow; label=label_low, marker=marker_low)
end

# bar_width is total bar width of each group, so bar_width=1 will pack all bar groups together
function errorbars_flipped!(data, low, high; bar_width=0.8, line=(color=:black, linewidth=6),
                    label_low="warm start", marker_low=(:hline, :black, 20),
                    label_high="upper bound", marker_high=(:hline, :dot, :black, 20))
    c, r = size(data)
    xgroup = repeat(0:r-1, inner=c)
    barindex = repeat(0:c-1, outer=r)
    width = bar_width / c
    xfirstbar = (1 - bar_width)/2
    x = xgroup .+ (xfirstbar + width/2 .+ width*barindex)
    ylow = low[:]
    yhigh = high[:]
    xlines = [x x fill(NaN, r*c)]'[:]
    ylines = [data[:] yhigh fill(NaN, r*c)]'[:]
    plot!(xlines, ylines; label="", line...)
    scatter!(x, yhigh, label=label_high, marker=marker_high)
    scatter!(x, ylow; label=label_low, marker=marker_low)
end

function barlines_flipped!(data; bar_width=0.8, label="recalculated to A", line=(color=:dodgerblue2, linewidth=2, linestyle=:dot))
    c, r = size(data)
    xgroup = repeat(0:r-1, inner=c)
    barindex = repeat(0:c-1, outer=r)
    width = bar_width / c
    xfirstbar = (1 - bar_width)/2
    x1 = xgroup .+ (xfirstbar .+ width*barindex)
    x2 = xgroup .+ (xfirstbar + width .+ width*barindex)
    xlines = [x1 x2 fill(NaN, r*c)]'[:]
    ylines = [data[:] data[:] fill(NaN, r*c)]'[:]
    plot!(xlines, ylines; label, line...)
end

#=
@recipe function f(groupby::StatsPlots.RecipesPipeline.GroupBy, args...; sortgroups=true)
    plt = plotattributes[:plot_object]
    group_length = maximum(union(groupby.group_indices...))
    if !(StatsPlots.RecipesPipeline.group_as_matrix(args[1]))
        for (i, glab) in enumerate(groupby.group_labels)
            @series begin
                label --> string(glab)
                idxfilter --> groupby.group_indices[i]
                for (key, val) in plotattributes
                    if StatsPlots.RecipesPipeline.splittable_attribute(plt, key, val, group_length)
                        :($key) := split_attribute(plt, key, val, groupby.group_indices[i])
                    end
                end
                args
            end
        end
    else
        g = args[1]
        if length(g.args) == 1
            x = zeros(Int, group_length)
            for indexes in groupby.group_indices
                x[indexes] = eachindex(indexes)
            end
            last_args = g.args
        else
            x = g.args[1]
            last_args = g.args[2:end]
        end
        x_u = sortgroups ? unique(sort(x)) : unique(x)
        x_ind = Dict(zip(x_u, eachindex(x_u)))
        for (key, val) in plotattributes
            if StatsPlots.RecipesPipeline.splittable_attribute(plt, key, val, group_length)
                :($key) := StatsPlots.RecipesPipeline.groupedvec2mat(x_ind, x, val, groupby)
            end
        end
        label --> reshape(groupby.group_labels, 1, :)
        typeof(g)((
            x_u,
            (StatsPlots.RecipesPipeline.groupedvec2mat(x_ind, x, arg, groupby, NaN) for arg in last_args)...,
        ))
    end
end
=#
