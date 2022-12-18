using MathOptInterface

include("modelBinPacking.jl")
include("heuristics.jl")

"""
    BnB1(solver, instance, heuristic)

    Branch and Bound (BnB) method entry point. Start BnB method to solve `instance`.
    Use `heuristic` to compute an upper bound of `m` and `solver` to solve each subproblems.

    # Arguments
    - `solver`: the solver Optimizer function.
    - `instance`: the bin packing instance.
    - `heuristic`: the heuristic used to compute an upper bound of `m` the number of bins.
    - `verbose` : prints stuff if enabled.
"""
function BnB1(solver, instance, heuristic, verbose::Bool = false)
    # Compute upper bound of m with heuristic procedure.
    # Also get indices of instance.w sorted by decreasing value (ind) and
    # sum of elements in instance.w (s)
    n::Int64 = length(instance.w)
    @time begin
        m::Int64, s::Int64, ind::Vector{Int64} = heuristic(instance)

        # Vector of residual capacities of the bins
        cbar::Vector{Int64} = fill(instance.C, m)

        # Nf[i] true if x[i] is free, false if x[i] is fixed
        Nf::Vector{Bool} = fill(true, m*n)

        # Initial problem : P0 is a tuple (model, x, y)
        P0 = modelBinPacking(solver, instance)

        # Results
        R0 = recBnB(P0, instance, s, m, m, 1, ind, cbar, Nf, 0, verbose)

        if R0[1] != -1 && R0[1] < m
            m = R0[1]
        end

        println("\n  Meilleure valeur obtenue : z = ", m, "\n")
    end

    println()
    return m, R0[2], R0[3]
end

function recBnB(P,
                instance,
                s::Int64,
                m::Int64,
                zbar::Int64,
                idxItem::Int64,
                ind::Vector{Int64},
                cbar::Vector{Int64},
                Nf::Vector{Bool},
                depth::Int64 = 0,
                verbose::Bool = false,
                bidx::Int64 = -1,
                bval::Int64 = -1)
    # Get the number of items
    n::Int64 = length(instance.w)
    # Solve the relaxation
    optimize!(P[1])

    if(verbose) println("----------- ", depth, " -----------") end
    if(verbose && depth > 0) println("x[", bidx, "] = ", bval) end
    if has_values(P[1])
        z::Float64 = objective_value(P[1])
        isint::Bool, k::Int64 = true, 1
        while k < length(P[1][:x]) && isint
            isint &= isinteger(value(P[1][:x][k]))
            k += 1
        end

        if(verbose) println("z = ", z) end
        if isint
            if(verbose) println("Noeud sondé (optimalité)") end
            return z, P[1][:x], P[1][:y]
        else
            L1::Float64 = s/instance.C

            if(verbose && depth > 0) println("L1 = ", L1) end
            if(verbose) println("Solution non réalisable") end

            if z > zbar
                if(verbose) println("Noeud sondé (dominance)") end
                return -1, nothing, nothing
            elseif L1 > z && depth > 0
                if(verbose) println("Noeud sondé (L1)") end
                return -1, nothing, nothing
            else
                # Choose item index and bin index
                i::Int64, item::Int64, bin::Int64, idx::Int64 = idxItem, -1, -1, -1
                while i <= n && item == -1 && idxItem <= n
                    j::Int64 = 1
                    while j <= m && bin == -1 && Nf[m * (i-1) + j]
                        if cbar[j] - instance.w[ind[i]] >= 0
                            # item, bin and branching variable indices
                            item, bin, idx = i, j, m * (i-1) + j
                            Nf[idx] = false
                        end
                        j += 1
                    end
                    i += 1
                end

                if depth > 0 && idxItem <= n
                    s -= instance.w[ind[idxItem]]
                end

                xbest, ybest = nothing, nothing
                if item == -1 || bin == -1 || idxItem > n
                    return -1, xbest, ybest
                end

                for v=0:1
                    # Subproblems
                    coupe = @constraint(P[1], P[1][:x][idx] == v)
                    R = recBnB(P, instance, s, m, zbar, idxItem+1, ind, cbar, Nf,
                               depth+1, verbose, idx, v)
                    delete(P[1], coupe)
                    if R[1] != -1 && R[1] < zbar
                        zbar, xbest, ybest = R
                    end
                    if v == 0
                        cbar[bin] -= instance.w[ind[item]]
                    end
                end

                Nf[idx] = true
                cbar[bin] += instance.w[ind[item]]
                return zbar, xbest, ybest
            end
        end
    else
        if(verbose) println("Noeud sondé (infaisabilité)") end
        return -1, nothing, nothing
    end
end
