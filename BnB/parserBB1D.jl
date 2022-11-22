# ==============================================================================
# parserBB1D.jl
# X. Gandibleux - Novembre 2022

# ------------------------------------------------------------------------------
#= Retourne les données d'une instance de BB1D au format OR-library
   The format of these data files is as follows:
     Number of test problems (P)
     For each test problem (p=1,...,P) in turn:
        Problem identifier
        Bin capacity, Number of items (n), best known solution
        For each item i (i=1,...,n): size of the item
=#

function loadBB1D(nomFichier::String)

    f = open(nomFichier,"r")

    # nombre d'instances dans le fichier
    P = parse(Int64, readline(f))
	# Reservation mémoire pour memoriser toutes les instances
	allInstances = Vector{instanceBB1D}(undef,P)

	# lecture de chaque instance
	for p in 1:P

		# identificateur de l'instance
		id = readline(f)
		# capacite du bin, nbr d'objets, best sol (info non conservee)
		C, n, zBest = parse.(Int64,split(readline(f)))
		# allocation memoire de l'instance
		allInstances[p] = instanceBB1D(id,C,Vector{Int64}(undef,n))

		# Lecture des objets; 1 ligne -> espace requis par l'objet
		for i in 1:n
			allInstances[p].w[i] = parse(Int64, readline(f))
		end
    end

    close(f)
    return allInstances
end
