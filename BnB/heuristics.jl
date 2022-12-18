"""
    NFD(instance)

    Compute an upper bound of `m` (the number of bins) for the given
    `instance` using the next-fit decreasing heuristic.

    # Arguments
    - `instance`: the bin packing instance.
"""
function NFD(instance)
    j::Int64, s::Int64 = 1, 0
    C::Int64 = instance.C
    n::Int64 = length(instance.w)
    ind = sortperm(instance.w, rev=true)
    bins::Vector{Set{Int64}} = [Set{Int64}()]

    for i=1:n
        s += instance.w[i]
        if C - instance.w[ind[i]] >= 0
            push!(bins[j], ind[i])
            C -= instance.w[ind[i]]
        else
            j += 1
            C = instance.C - instance.w[ind[i]]
            push!(bins, Set{Int64}())
            push!(bins[j], ind[i])
        end
    end

    return length(bins), s, ind
end

"""
    FFD(instance)

    Compute an upper bound of `m` (the number of bins) for the given
    `instance` using the first-fit decreasing heuristic.

    # Arguments
    - `instance`: the bin packing instance.
"""
function FFD(instance)
    j::Int64, s::Int64 = 1, 0
    n::Int64 = length(instance.w)
    ind = sortperm(instance.w, rev=true)
    capa::Vector{Int64} = [instance.C]
    bins::Vector{Set{Int64}} = [Set{Int64}()]

    for i=1:n
        s += instance.w[i]
        stored::Bool = false
        # Trouver un bin avec assez de place
        k = 1
        while k <= j && !stored
            if capa[k] - instance.w[ind[i]] >= 0
                push!(bins[k], ind[i])
                capa[k] -= instance.w[ind[i]]
                stored = true
            end
            k += 1
        end

        # Si aucun bin n'a été trouvé, on met l'objet dans un nouveau bin
        if !stored
            j += 1
            push!(bins, Set{Int64}())
            push!(bins[j], ind[i])
            push!(capa, instance.C - instance.w[ind[i]])
        end
    end

    return length(bins), s, ind
end

"""
    BFD(instance)

    Compute an upper bound of `m` (the number of bins) for the given
    `instance` using the best-fit decreasing heuristic.

    # Arguments
    - `instance`: the bin packing instance.
"""
function BFD(instance)
    j::Int64, s::Int64 = 1, 0
    n::Int64 = length(instance.w)
    ind = sortperm(instance.w, rev=true)
    capa::Vector{Int64} = [instance.C]
    bins::Vector{Set{Int64}} = [Set{Int64}()]

    for i=1:n
        s += instance.w[i]
        kmin, capamin = -1, typemax(Int64)
        # Trouver un bin avec assez de place
        for k=1:j
            if capa[k] < capamin && capa[k] >= instance.w[ind[i]]
                kmin, capamin = k, capa[k]
            end
        end

        # Si le bin de plus petite capacité résiduelle n'a pas assez de place
        # => nouveau bin
        if capamin < instance.w[ind[i]] || kmin == -1
            j += 1
            push!(bins, Set{Int64}())
            push!(bins[j], ind[i])
            push!(capa, instance.C - instance.w[ind[i]])
        else
            push!(bins[kmin], ind[i])
            capa[kmin] -= instance.w[ind[i]]
        end
    end

    return length(bins), s, ind
end
