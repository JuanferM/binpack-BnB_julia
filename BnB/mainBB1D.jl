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
        z, x, y = BnB1(GLPK.Optimizer, instance, BFD)
        if x !== nothing
            println(value.(x))
            println(length(x))
            println(value.(y))
        end

        # binpacking, x, y = modelBinPacking(GLPK.Optimizer, instance, 5, true, true)
        #
        # println("\nOptimisation...")
        # optimize!(binpacking)
        # println("\nRésultats")
        # println(solution_summary(binpacking))
        # println(value.(x))
        # println(length(x))
        # println(value.(y))
        #
        # println(node_count(binpacking))
    end

    return nothing
end

main("didactic.txt")
