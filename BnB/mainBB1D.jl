# ==============================================================================
# mainBB1D.jl
# X. Gandibleux - Novembre 2022
# Modifié par Juanfer MERCIER & Tobias SANVOISIN - Novembre 2022

using JuMP
using GLPK

include("datastrucBB1D.jl")
include("parserBB1D.jl")
include("modelBinPacking.jl")
include("heuristics.jl")

function main(fname::String, use_solver::Bool = false)
    data::Vector{instanceBB1D} = loadBB1D("instances/" * fname)

    for instance in data
        println(instance.id)
        if use_solver
            binpacking, x, y = modelBinPacking(GLPK.Optimizer, instance, true)
            println("\nOptimisation...")
            optimize!(binpacking)
            println("\nRésultats")
            println(solution_summary(binpacking))
        else
            print(NFD(instance))
        end
    end

    return nothing
end

main("didactic.txt")
