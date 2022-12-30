# ==============================================================================
# mainBB1D.jl
# X. Gandibleux - Novembre 2022
# Modifié par Juanfer MERCIER & Tobias SANVOISIN - Novembre 2022

using JuMP
using GLPK
using Gurobi
using MathOptInterface

include("datastrucBB1D.jl")
include("parserBB1D.jl")
include("modelBinPacking.jl")
include("heuristics.jl")
include("BnB.jl")

"""
    main(fname)

    Main entry point of the method, run method on the instances stored in `fname`.

    # Arguments
    - `fname`: the filename of the file containing the instances.
"""
function main(fname::String)
    data::Vector{instanceBB1D} = loadBB1D("instances/" * fname)

    for instance in data
        println(instance.id)
        # Use BnB with GLPK on instance with BFD heuristic.
        # Verbose disabled. Timeout at 45 seconds.
        # By default, we use L2 and minimal cover cuts.
        # z, x, y = BnB(GLPK.Optimizer, instance, BFD, false, 45.0)
        z, x, y = BnB(Gurobi.Optimizer, instance, BFD, false, 45.0)

        # m::Int64, _, _ = FFD(instance)
        # binpacking, x, y = modelBinPacking(GLPK.Optimizer, instance, m, false, false)
        # set_time_limit_sec(binpacking, 45.0)
        #
        # println("\nOptimisation...")
        # optimize!(binpacking)
        # println("\nRésultats")
        # println(solution_summary(binpacking))
    end

    return nothing
end

main("N1C1W1.txt")
