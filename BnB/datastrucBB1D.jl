# ==============================================================================
# datastrucBB1D.jl
# X. Gandibleux - Novembre 2022

# ------------------------------------------------------------------------------
# Structures contenant les données d'un problème de bin-packing mono-dimensionnel
struct instanceBB1D
	id :: String        # identifiant de l'instance
	C  :: Int64         # capacite d'un bin
	w  :: Vector{Int64} # vecteur de l'espace requis par les n objet
end
