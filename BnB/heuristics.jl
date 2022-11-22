# next-fit decreasing
function NFD(instance)
    n::Int64 = length(instance.w)
    bins::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef, n)
    for i=1:n
        bins[i] = []
    end
end

# first-fit decreasing
function FFD(instance)

end

# best-fit decreasing
function BFD(instance)

end
