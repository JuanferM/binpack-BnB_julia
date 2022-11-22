# ==============================================================================
# mainBB1D.jl
# X. Gandibleux - Novembre 2022
# Modifié par Juanfer MERCIER & Tobias SANVOISIN - Novembre 2022

using JuMP
using GLPK

include("datastrucBB1D.jl")
include("parserBB1D.jl")
include("modelBinPacking.jl")

function main(fname::String)
    data::Vector{instanceBB1D} = loadBB1D("instances/" * fname)

    for instance in data
        println(instance.id)
        binpacking, x, y = modelBinPacking(GLPK.Optimizer, instance, true)
        println("\nOptimisation...")
        optimize!(binpacking)
        println("\nRésultats")
        println(solution_summary(binpacking))
    end

    return nothing
end

main("didactic.txt")
