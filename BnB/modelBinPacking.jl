using JuMP

"""
    modelBinPacking(solver, instance[, bound, P1, binary])

    Create the bin packing model corresponding to _instance_ using
    `solver` as the optimizer. When `P1` is true we model against problem P1 (see
    report for further details), when `P1` is false we model against problem P2.
    If `binary` is true we consider boolean variables, if `binary` is false we consider
    free variables.

    # Arguments
    - `solver`: the solver Optimizer function.
    - `instance`: the bin packing instance.
    - **optional** `bound`: upper bound of m. By default, m = n.
    - **optional** `P1`: true if modeling with P1 (P2 otherwise and by default).
    - **optional** `binary`: true if modeling with binary variables (free variables otherwise).
"""
function modelBinPacking(solver, instance, bound::Int64 = -1, P1::Bool = false, binary::Bool = false)
    n::Int64 = length(instance.w)
    m::Int64 = (bound == -1) ? n : bound

    # Création du model
    model = direct_model((solver)())
    if solver == Gurobi.Optimizer
    	set_optimizer_attribute(model, "OutputFlag", 0)
    end

    # Définition des variables
    @variable(model, 0 <= x[1:(n*m)] <= 1, binary=binary)
    @variable(model, 0 <= y[1:m] <= 1, binary=binary)

    # Définition de l'objectif (la somme des yj à minimiser)
    @objective(model, Min, sum(y[j] for j=1:m))

    # Définition des contraintes
    @constraint(model, assign[i=1:n], sum(x[(i-1)*m+j] for j=1:m) == 1)
    if(P1)
        @constraint(model, capa[j=1:m], sum(x[(i-1)*m+j]*instance.w[i] for i=1:n) <= instance.C * y[j])
    else
        # Rajout de la contrainte additionnelle
        @constraint(model, capa[j=1:m], sum(x[(i-1)*m+j]*instance.w[i] for i=1:n) <= instance.C)
        @constraint(model, ranger[i=1:n, j=1:m], x[(i-1)*m+j] <= y[j])
    end

    return model, x, y
end
