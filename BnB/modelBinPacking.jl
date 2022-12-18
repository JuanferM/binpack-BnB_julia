using JuMP

"""
    modelBinPacking(solver, instance[, P1, binary])

    Create the bin packing model corresponding to _instance_ using
    `solver` as the optimizer. When `P1` is true we model against problem P1 (see
    report for further details), when `P1` is false we model against problem P2.
    If `binary` is true we consider boolean variables, if `binary` is false we consider
    free variables.

    # Arguments
    - `solver`: the solver Optimizer function.
    - `instance`: the bin packing instance.
    - **optional** `P1`: true if modeling with P1 (P2 otherwise).
    - **optional** `binary`: true if modeling with binary variables (free variables otherwise).
    - **optional** `bound`: upper bound of m. By default, m = n.
"""
function modelBinPacking(solver, instance, P1::Bool = true, binary::Bool = false, bound = -1)
    n = length(instance.w)
    m = (bound == -1) ? n : bound

    # Création du model
    model = direct_model((solver)())

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
