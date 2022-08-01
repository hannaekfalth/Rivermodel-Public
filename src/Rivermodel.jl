module Rivermodel

using JuMP, Gurobi, Ipopt, AxisArrays, Plots, UnPack, FileIO, Statistics, Polynomials

export runmodel, read_inputdata, read_raw_data, exploredata, printmetrics, compareruns

struct Upstream
    name::Symbol
    dischargedelay::Int     # hours
    spilldelay::Int         # hours
    flowshare::Float64      # fraction of upstream flow that reaches plant
                            # (< 1 for the first plant in a branch after a river fork) 
end
Upstream(name, ddelay, sdelay) = Upstream(name, ddelay, sdelay, 1.0)

struct Plant
    name::Symbol
    capacity::Float64           # MW
    reservoir::Float64          # HE
    reportedhead::Float64       # m
    reservoirhigh::Float64      # m
    reservoirlow::Float64       # m
    upstream::Vector{Upstream}  # vector of upstream plants
                                # (more than one if two river branches converge upstream)
end

Plant(name, cap, res, head, rh, rl, p::Upstream) =
        Plant(name, cap, res, head, rh, rl, [p])
Plant(name, cap, res, head, rh, rl, upstream::Upstream...) =
        Plant(name, cap, res, head, rh, rl, [upstream...])

AxisArrays.AxisArray(f::Float64, axes...) = AxisArray(fill(f, length.(axes)), axes...)
macro aa(e)
    esc(:(AxisArray($e...)))
end

const DATAFOLDER = Sys.iswindows() ? "C:/Skelleftedata" : "/mnt/c/Skelleftedata"

include("rivernetworks.jl")
include("skellefteinput.jl")
include("helpfunctions.jl")
include("jumpmodel.jl")
include("output.jl")
include("plots.jl")
include("rawdata.jl")
include("exploredata.jl")
include("compareruns.jl")
include("aggmodel.jl")

# USAGE:
# type = [:LP, :NLP, :MIP]
# power = ["E constant head", "E taylor", "bilinear HeadE", "aggregated"]
# e = ["cv segments origo", "cv segments noseg", "constant eta", "ncv segments rampseg", "ncv segments zeroseg", "cv poly noseg", "ncv poly rampseg", "ncv poly zeroseg"]
# head = [:regression, :constantmean, :constantmax, :reservoirdiff, :singlereservoir, :aggregated]
# start = (type=:LP, power="E taylor", e="convex segments origo") or just (type=:LP,) to use same power & e as main run (but don't forget the last comma).
# recalc [optional], for recalculating an LP run when there is no run 2. default: recalc = (type=:NLP, power="bilinear HeadE", e="ncv poly rampseg") 
# regression variable order: constant, discharge, discharge lags (1,2,3), forebay, downstream_forebay

function runmodel(model_year; weeks=1:52, env_con=true, type=:LP, power="taylor E", e="cv segments origo", head=:standardbasic,
            start::NamedTuple=(;), recalc::NamedTuple=(;), save_variables=false, runfamily = "", regressionvars=Bool[1,1,1,1,1,1,1], silent=true)

    if head == :standardbasic
        regressionvars = Bool[1,1,0,0,0,0,0]
    elseif head == :standardplus
        regressionvars = Bool[1,1,0,0,0,0,1]
    elseif head == :standardall
        regressionvars = Bool[1,1,1,1,1,1,1]
    elseif !(head in [:regression, :constantmean, :constantmax, :reservoirdiff, :singlereservoir])
        error("$head is not a recognized head type.")
    end

    isaggregated = (power == "aggregated")
    start = isempty(start) ? start : (; power, e, head, start...) # use main power & e arguments as defaults (so no need to repeat them if identical)
    run2args = (; type, power, e, head) 
    run1args = isempty(start) ? run2args : (type=start.type, power=start.power, e=start.e, head=start.head)
    recalcargs = isempty(start) ? (type=:NLP, power="bilinear HeadE", e="ncv poly rampseg", head, recalc...) : run2args

    @time params = read_inputdata(model_year, head; weeks, env_con, regressionvars, silent)
    @time results = isaggregated ? buildmodel_agg(params, model_year) : buildmodel(params, model_year; run1args...)

    rivermodel = results.rivermodel

    println("Solving model...")
    
    firsttype = isempty(start) ? type : start.type
    setsolver(rivermodel, (firsttype == :NLP) ? :ipopt : :gurobi)
    get(start, :type, :novalue) == :MIP && set_optimizer_attribute(rivermodel, "SolutionLimit", 1)
    optimize!(rivermodel)

    status = termination_status(rivermodel)
    status != MOI.OPTIMAL && @warn "The solver did not report an optimal solution."
    println("\nSolve status: $status")

    printbasicresults(params, results; isaggregated, recalcargs..., recalculate=true)
    runtype = (type == :LP || isempty(start)) ? "main" : "start"
    save_variables && savevariables(model_year, results, runtype, runfamily, solve_time(rivermodel); run1args...)

    if type == :LP || isempty(start)
        return rivermodel, params, results
    end

    
    println("\n\nBuilding second model (because modifying JuMP models is super slow)...")
    @time results2 = buildmodel(params, model_year; run2args...)
    rivermodel2 = results2.rivermodel

    println("\nSetting variable start values to LP result...")
    vars = filter(x -> name(x)[1] != 'Z', all_variables(rivermodel)) # if model 1 is a MIP, ignore the binary variables
    vars2 = all_variables(rivermodel2)
    set_start_value.(vars2, value.(vars))
    set_start_values!(params, results, results2; run2args...)
                    
    println("\nSolving model with start values...")
    setsolver(rivermodel2, (type == :NLP) ? :ipopt : :gurobi)
    optimize!(rivermodel2)

    status = termination_status(rivermodel2)
    if type == :LP && status != MOI.OPTIMAL
        @warn "The solver did not report an optimal solution."
    end
    println("\nSolve status: $status")

    printbasicresults(params, results2; run2args..., recalculate=false, isaggregated)
    save_variables && savevariables(model_year, results2, "main", runfamily, solve_time(rivermodel2); run2args..., start)

    return rivermodel2, params, results2
end


function setsolver(model, solver)
    if solver == :gurobi
        # Workaround for bug:
        # https://discourse.julialang.org/t/how-can-i-clear-solver-attributes-in-jump/57953
        optimizer = optimizer_with_attributes(Gurobi.Optimizer, "Threads" => 4, "Method" => 2, "Presolve" => 2, "PreSparsify" => 1, "Cuts" => 2, 
            "nonconvex" => 2, "crossover" => 0, "MIPGap" => 5e-6, "DisplayInterval" => 1)
        set_optimizer(model, optimizer)
                # "GURO_PAR_DUMP" => 1  # output an MPS file for external solving
                # "Heuristics" => 0.9,  # fraction of total MIP runtime devoted to heuristics 
                # "PreMIQCPForm" => 2,  # -1, 0, 1, 2
                # "MIPFocus" => 3,      # 1 finding feasible, 2 proving optimality, 3 move bound
                # "Cuts" => 1,          # -1 auto, 0 off, 1 on, 2 aggressive, 3 very aggressive
                # "Presolve" => 1,      # -1 auto, 0 off, 1 on, 2 aggressive
                # "PreSparsify" => 1,   # -1 auto, 0 off, 1 on
                # "SolutionLimit" => 1,
                # "NoRelHeurTime" => 3600*24*3,
    elseif solver == :ipopt
        optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 100000, "mu_strategy" => "adaptive")
                                            #"mumps_mem_percent" => 500, #1000 is default
                                            #"mumps_pivtol" => 1e-2, # range 0 to 1, 10^-6 is default
                                            #"tol" => 5e-7
        set_optimizer(model, optimizer)
        Sys.islinux() && set_optimizer_attributes(model, "linear_solver" => "ma86", "ma86_scaling" => "mc64")  # none, mc64 (default), mc77
        # Uncomment the following line to make IPOPT not go infeasible on iteration 1 (the first two options are enough). This is only for testing though:
        # Since IPOPT is an interior point solver warm starting is not a good idea (it needs to go into the interior and not remain on the bounds).
        # set_optimizer_attributes(model, "warm_start_init_point" => "yes", "warm_start_bound_push" => 1e-9, "warm_start_bound_frac" => 1e-9,
        #         "warm_start_slack_bound_frac" => 1e-9, "warm_start_slack_bound_push" => 1e-9, "warm_start_mult_bound_push" => 1e-9)
    else
        @error "No solver named $solver."
    end
end

end #module
