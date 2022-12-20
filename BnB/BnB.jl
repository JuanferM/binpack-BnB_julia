using MathOptInterface

include("modelBinPacking.jl")
include("heuristics.jl")

mutable struct Data
    z::Float64
    node_count::Int64
    solution_found::Bool
    solution_count::Int64
    start_time::Float64
end

"""
    BnB1(solver, instance, heuristic)

    Branch and Bound (BnB) method entry point. Start BnB method to solve `instance`.
    Use `heuristic` to compute an upper bound of `m` and `solver` to solve each subproblems.

    # Arguments
    - `solver`: the solver Optimizer function.
    - `instance`: the bin packing instance.
    - `heuristic`: the heuristic used to compute an upper bound of `m` the number of bins.
    - **optional** `verbose` : prints stuff if enabled.
    - **optional** `timeout` : search is stopped if timeout time is exceeded.
"""
function BnB1(solver,
              instance,
              heuristic,
              verbose::Bool = false,
              timeout::Float64 = typemax(Float64))
    # Compute upper bound of m with heuristic procedure.
    # Also get indices of instance.w sorted by decreasing value (ind) and
    # sum of elements in instance.w (s)
    n::Int64 = length(instance.w)
    @time begin
        m::Int64, s::Int64, ind::Vector{Int64} = heuristic(instance)

        # Vector of residual capacities of the bins
        cbar::Vector{Int64} = fill(instance.C, m)

        # Best z value and counter on the number of node explored
        data::Data = Data(m+0.0, 0, false, 0, time())

        # Nf[i] true if x[i] is free, false if x[i] is fixed
        Nf::Vector{Bool} = fill(true, m*n)

        # Initial problem : P0 is a tuple (model, x, y)
        P0 = modelBinPacking(solver, instance, m)

        # Results
        R0 = recBnB(P0, instance, s, m, 0, 0, -1, data, ind,
                    cbar, Nf, 0, verbose, false, timeout)

        println("\n  Nombre de noeuds explorés : ", data.node_count)
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
                fallback::Bool = false,
                timeout::Float64 = typemax(Float64))
    L1::Int64, idx::Int64, n::Int64 = -1, -1, length(instance.w)
    coupe::Union{ConstraintRef, Nothing}, modified::Bool = nothing, false

    println(time()-data.start_time)

    # Reached a leaf or exceeded time so leave.
    if(item < 0 || item > n || bin < 0 || bin > m
       || depth > n*m+1 || time()-data.start_time > timeout)
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

        L1 = Int64(ceil(s/instance.C)) # Compute L1
        idx = m * (item-1) + bin # Compute branching variable index

        # Add constraints
        Nf[idx] = false
        coupe = @constraint(P[1], P[1][:x][idx] == val)

        if(verbose)
            println("x[$(idx)] = $(val)")
            println("L1 = $(L1)")
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
            elseif floor(z) < L1 && data.solution_found # Worst than L1 bound
                if(verbose) println("Noeud sondé (L1)") end
            else
                # compute branching index and branch
                i::Int64, j::Int64, mincapa::Int64, capa::Int64 = item+1, 1, typemax(Int64), -1
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
                                  ind, cbar, Nf, depth+1, verbose, fallback, timeout)
                    if subR[1] != -1
                        if(subR[1] < Rz || Rz == -1) Rz, Rx, Ry = subR end
                    else
                        fallback = true
                    end
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
        Nf[idx] = true
        if(coupe !== nothing) delete(P[1], coupe) end
    end

    # Return results
    return Rz, Rx, Ry
end
