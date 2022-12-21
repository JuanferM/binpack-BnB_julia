# ==============================================================================
# mainBB1D.jl
# X. Gandibleux - Novembre 2022
# Modifié par Juanfer MERCIER & Tobias SANVOISIN - Novembre 2022

using JuMP
using GLPK
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
        z, x, y = BnB1(GLPK.Optimizer, instance, BFD, false, 45.0, false)

        # m::Int64, _, _ = BFD(instance)
        # binpacking, x, y = modelBinPacking(GLPK.Optimizer, instance, m, false, true)
        # set_time_limit_sec(binpacking, 45.0)
        #
        # println("\nOptimisation...")
        # optimize!(binpacking)
        # println("\nRésultats")
        # println(solution_summary(binpacking))
        # println(node_count(binpacking))
    end

    return nothing
end

main("N1C1W1.txt")
