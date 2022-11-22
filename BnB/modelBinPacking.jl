using JuMP

# Paramètres
#   solveur  : le solveur choisie pour traiter le problème
#   instance : une instance du bin packing
#   P1       : booléen à true si on choisit de résoudre le problème P1
#       si P1 est à false alors on résout le problème P2
#   bin      : booléen à true si on résout le problème en variables binaires
# Valeur de retour
#   un modèle JuMP du binpacking presque complet
function modelBinPacking(solveur, instance, bin::Bool = true, P1::Bool = true)
    n = length(instance.w)
    m = n

    # Création du model
    binpacking = Model()

    # Attribution du solveur
    if !isnothing(solveur)
        set_optimizer(binpacking, solveur)
    end

    # Définition des variables
    @variable(binpacking, 0 <= x[1:n, 1:m] <= 1, binary=bin)
    @variable(binpacking, 0 <= y[1:m] <= 1, binary=bin)

    # Définition de l'objectif (la somme des yj à minimiser)
    @objective(binpacking, Min, sum(y[j] for j=1:m))

    # Définition des contraintes
    @constraint(binpacking, assign[i=1:n], sum(x[i,j] for j=1:m) == 1)
    if(P1)
        @constraint(binpacking, capa[j=1:m], sum(x[i,j]*instance.w[i] for i=1:n) <= instance.C * y[j])
    else
        # Rajout de la contrainte additionnelle
        @constraint(binpacking, capa[j=1:m], sum(x[i,j]*instance.w[i] for i=1:n) <= instance.C)
        @constraint(binpacking, ranger[i=1:n, j=1:m], x[i, j] <= y[j])
    end

    return binpacking, x, y
end
