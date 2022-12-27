using MathOptInterface

include("modelBinPacking.jl")
include("heuristics.jl")

mutable struct Data
    z::Float64
    node_count::Int64
    solution_found::Bool
    solution_count::Int64
    global_L2::Int64
    start_time::Float64
    timed_out::Bool
end

"""
    BnB1(solver, instance, heuristic[, use_L2, only_free, verbose, timeout])

    Branch and Bound (BnB) method entry point. Start BnB method to solve `instance`.
    Use `heuristic` to compute an upper bound of `m` and `solver` to solve each subproblems.

    # Arguments
    - `solver`: the solver Optimizer function.
    - `instance`: the bin packing instance.
    - `heuristic`: the heuristic used to compute an upper bound of `m` the number of bins.
    - **optional** `verbose` : prints stuff if enabled.
    - **optional** `timeout` : search is stopped if timeout time is exceeded.
    - **optional** `use_L2` : use L2 bound instead of L1.
    - **optional** `only_free` : only use the free variables weights to compute L2.
    - **optional** `use_cut` : use minimal cover cuts
"""
function BnB(solver,
             instance,
             heuristic,
             verbose::Bool = false,
             timeout::Float64 = typemax(Float64),
             use_L2::Bool = true,
             only_free::Bool = false,
             use_cut::Bool = true)
    # Compute upper bound of m with heuristic procedure.
    # Also get indices of instance.w sorted by decreasing value (ind) and
    # sum of elements in instance.w (s)
    n::Int64 = length(instance.w)
    @time begin
        m::Int64, s::Int64, ind::Vector{Int64} = heuristic(instance)

        # Vector of residual capacities of the bins
        cbar::Vector{Int64} = fill(instance.C, m)

        # Best z value and counter on the number of node explored
        data::Data = Data(m+0.0, 0, false, 0, -1, time(), false)

        # Nf[i] true if x[i] is free, false if x[i] is fixed
        Nf::Vector{Bool} = fill(true, m*n)

        # Initial problem : P0 is a tuple (model, x, y)
        P0 = modelBinPacking(solver, instance, m)

        # Compute L2 if necessary (only once if only_free is false)
        # If only_free is true then L2 is computed at every node and
        # global_L2 is not used
        if use_L2 && !only_free
            data.global_L2 = getL2(instance, m, ind, Nf)
        end

        # Results
        R0 = recBnB(P0, instance, s, m, 0, 0, -1, data, ind, cbar,
                    Nf, 0, verbose, timeout, use_L2, only_free, use_cut)
        println("\n")
        if(data.timed_out) println("  Expiration du temps accordé!") end
        println("  Nombre de noeuds explorés : ", data.node_count)
        println("  Nombre de solutions réalisables trouvées : ", data.solution_count)
        println("  Meilleure valeur obtenue : z = ", data.z, "\n")
    end

    println()
    return R0
end

"""
    checkFeasibility(x, y)

    Check if x and y are vectors of integers (if every element of x and y is either equal to 0
    or 1). x and y should be vectors of variable references (VariableRef). Because of that, we
    also create two vectors of floats while doing the checking that we return if x and y are true.

    # Arguments
    - `x`: the Vector{VariableRef} representing x.
    - `y`: the Vector{VariableRef} representing y.
"""
function checkFeasibility(x, y)
    if(x === nothing || y === nothing) return false, x, y end
    lx::Int64, ly::Int64 = length(x), length(y)
    nx::Vector{Float64}, ny::Vector{Float64} = [], []
    isintx::Bool, isinty::Bool, k::Int64, l::Int64 = true, true, 1, 1
    while k <= lx && isintx && isinty
        isintx &= isinteger(value(x[k]))
        push!(nx, value(x[k]))
        if l <= ly
            isinty &= isinteger(value(y[l]))
            push!(ny, value(y[l]))
        end
        k += 1
        l += 1
    end

    if(isintx && isinty) return true, nx, ny
    else return false, nothing, nothing end
end

"""
    getL2(instance, m, ind, Nf[, only_free])

    Compute L2 bound.

    # Arguments
    - `instance`: the bin packing instance.
    - `m`: bound on m obtained via heuristic algorithm.
    - `ind`: indices of the items ordered by decreasing weight.
    - `Nf`: boolean vector. An element of Nf is true if the variable of corresponding index is free.
    - **optional** `only_free`: take only free variables into account.
"""
function getL2(instance,
               m::Int64,
               ind::Vector{Int64},
               Nf::Vector{Bool},
               only_free::Bool = false)
    n::Int64 = length(instance.w)
    L2::Int64, C::Int64 = -1, instance.C

    for α::Int64 = 0:C/2
        lJ1::Int64, lJ2::Int64, sJ2::Int64, sJ3::Int64 = 0, 0, 0, 0

        for i=1:n
            compute::Bool, j::Int64 = true, (i-1)*m+1
            while j <= ((i-1)*m+m) && compute && only_free
                compute = Nf[j]
                j += 1
            end

            if compute
                if instance.w[ind[i]] > C - α
                    lJ1 += 1
                elseif C - α >= instance.w[ind[i]] > C/2
                    lJ2, sJ2 = lJ2+1, sJ2+instance.w[ind[i]]
                elseif C/2 >= ind[i] >= α
                    sJ3 += instance.w[ind[i]]
                end
            end
        end

        L2 = max(L2, lJ1 + lJ2 + max(0, ceil((sJ3 - (lJ2*C - sJ2))/C)))
    end

    return L2
end

"""
    recBnB(P, instance, s, m, item, bin, val, data, ind, cbar, Nf[, depth, verbose, fallback, timeout])

    Branch and Bound (BnB) branching method. Branching is done via recursive calls using the data at
    hand.

    # Arguments
    - `P`: the base problem that we'll be adding constraints to.
    - `instance`: the bin packing instance.
    - `s`: the sum of the items' weights that aren't in a bin.
    - `m`: bound on m obtained via heuristic algorithm.
    - `item`: index of the item to assign to the bin `bin`.
    - `bin`: index of the bin that `item` is assigned to.
    - `val`: the value that the variable (m * (item-1) + j) will be receiving.
    - `data`: holds data such as the value of the best solution, the number of nodes explored, etc.
    - `ind`: indices of the items ordered by decreasing weight.
    - `cbar`: residual capacities of the bins.
    - `Nf`: boolean vector. An element of Nf is true if the variable of corresponding index is free.
    - **optional** `depth`: depth of the current node.
    - **optional** `verbose`: prints stuff if enabled.
    - **optional** `timeout`: If BnB takes more than `timeout` seconds to finish, the search is stopped.
    - **optional** `use_L2` : use L2 bound instead of L1.
    - **optional** `only_free` : only use the free variables weights to compute L2.
    - **optional** `use_cut` : use minimal cover cuts
    - **optional** `fallback`: An item wasn't able to fit into a bin without infeasibility. `fallback`
        is true if when that happens so we select that same item to put it into another bin.
"""
function recBnB(P,
                instance,
                s::Int64,
                m::Int64,
                item::Int64,
                bin::Int64,
                val::Int64,
                data::Data,
                ind::Vector{Int64},
                cbar::Vector{Int64},
                Nf::Vector{Bool},
                depth::Int64 = 0,
                verbose::Bool = false,
                timeout::Float64 = typemax(Float64),
                use_L2::Bool = true,
                only_free::Bool = false,
                use_cut::Bool = true,
                fallback::Bool = false)
    L::Int64, index::Int64, n::Int64 = -1, -1, length(instance.w)
    cont::Union{ConstraintRef, Nothing}, modified::Bool = nothing, false

    # Reached a leaf or exceeded time so leave.
    if(time()-data.start_time > timeout)
        data.timed_out = true
    end
    if(item < 0 || item > n || bin < 0 || bin > m || depth > n*m+1 || data.timed_out)
        return -1, nothing, nothing
    end
    if(verbose) println("----------- $(depth) -----------") end

    # If not first node then modify s and cbar accordingly
    # Also add constraints
    if depth > 0
        modified = true
        if val == 1
            s -= instance.w[ind[item]]
            cbar[bin] -= instance.w[ind[item]]
        end

        # Compute bound (L1 or L2 according to use_L2 and only_free)
        if !use_L2
            L = Int64(ceil(s/instance.C))
        elseif use_L2 && !only_free
            L = data.global_L2
        else
            L = getL2(instance, m, ind, Nf, only_free)
        end
        index = m * (item-1) + bin # Compute branching variable index

        # Add constraints
        Nf[index] = false
        cont = @constraint(P[1], P[1][:x][index] == val)

        if(verbose)
            println("x[$(index)] = $(val)")
            println((use_L2 ? "L2" : "L1" ), " = $(L)")
            println("cbar = $(cbar)")
        end
    end

    optimize!(P[1]) # Compute the relaxation problem
    data.node_count += 1 # Increment node counter
    # Results
    Rz::Float64 = -1.0
    Rx::Union{Vector{Float64}, Nothing}, Ry::Union{Vector{Float64}, Nothing} = nothing, nothing

    if termination_status(P[1]) != INFEASIBLE
        z::Float64 = objective_value(P[1])
        isint::Bool, x, y = checkFeasibility(P[1][:x], P[1][:y])

        if(verbose) println("z = $(z)") end
        if isint
            if(verbose) println("Noeud sondé (optimalité)") end
            data.solution_count += 1
            if(z < data.z || !data.solution_found)
                data.z = z
                data.solution_found = true
                Rz, Rx, Ry = z, x, y
            end
        else
            if(verbose) println("Solution non réalisable") end

            if z > data.z # Worst than best solution until now
                if(verbose) println("Noeud sondé (dominance)") end
            elseif floor(z) < L && data.solution_found # Worst than the bound (L1 or L2)
                if(verbose) println("Noeud sondé (", (use_L2 ? "L2" : "L1"), ")") end
            else
                # compute cuts
                cuts::Vector{ConstraintRef} = []
                if use_cut
                    for j=1:m # compute minimal cover
                        i::Int64, s::Int64= 1, 0
                        cover::Vector{VariableRef}, cvi::Vector{Int64} = [], []
                        while i <= n && s <= instance.C
                            idx::Int64 = m * (ind[i]-1) + j
                            if !Nf[idx] # Item is in bin
                                s += instance.w[ind[i]]
                                push!(cover, P[1][:x][idx])
                                push!(cvi, ind[i])
                            end
                            i += 1
                        end

                        if length(cover) != 0 # if we have a cover
                            # add the cut to the problem
                            cut::ConstraintRef = @constraint(P[1], sum(cover) <= length(cover)-1)
                            push!(cuts, cut) # add the cut to the list
                        end
                    end
                end

                # compute branching index and branch
                i, j::Int64, mincapa::Int64, capa::Int64 = item+1, 1, typemax(Int64), -1
                if(fallback) i -= 1; fallback = false end
                newitem::Int64, newbin::Int64 = -1, -1
                while i <= n && newitem == -1
                    j = 1
                    while j <= m
                        capa = cbar[j] - instance.w[ind[i]]
                        if capa >= 0 && capa < mincapa && Nf[m * (i-1) + j]
                            newitem, newbin, mincapa = i, j, capa
                        end
                        j += 1
                    end
                    i += 1
                end

                # Branch and get best results
                for v=1:-1:0
                    subR = recBnB(P, instance, s, m, newitem, newbin, v, data,
                                  ind, cbar, Nf, depth+1, verbose, timeout,
                                  use_L2, only_free, fallback)
                    if subR[1] != -1
                        if(subR[1] < Rz || Rz == -1) Rz, Rx, Ry = subR end
                    else
                        fallback = true
                    end
                end

                # Remove cover cuts from the problem
                for cut in cuts
                    delete(P[1], cut)
                end
            end
        end
    else
        if(verbose) println("Noeud sondé (infaisabilité)") end
    end

    # If s and cbar has been modified get everything back to normal
    # Also remove constraints
    if modified
        if val == 1
            s += instance.w[ind[item]] # Doesn't change anything but whatever...
            cbar[bin] += instance.w[ind[item]] # Important
        end

        # remove constraints
        Nf[index] = true
        if(cont !== nothing) delete(P[1], cont) end
    end

    # Return results
    return Rz, Rx, Ry
end
