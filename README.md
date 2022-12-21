# binpacking_julia
Branch-and-bound algorithm and heuristics for the 1D bin packing problem in Julia.

The branch-and-bound algorithm (`BnB.jl`) follows the guidelines of Eilon and Christofides (1971) with
the addition of the L2 bound presented by Marthelo and Toth. The algorithm was tested on one instance 
of the bin packing problem considered in E. Falkenauer (1994) "A Hybrid Grouping Genetic Algorithm for
Bin Packing". This data file is binpack1 and was contributed by E. Falkenauer.

The format of this data file is as follows:
Number of test problems (P)
For each test problem (p=1,...,P) in turn:
   Problem identifier
   Bin capacity, Number of items (n), Number of bins in the current
                                      best known solution
   For each item i (i=1,...,n): size of the item
   
There also instances  used by A. Scholl, R. Klein, and C. Jürgens in Bison: a fast hybrid procedure 
for exactly solving the one-dimensional bin packing problem. Those instances were compiled into data
files in the same format as the data files considered by E. Falkenauer. Basically, all instances with
the same characteristics were put into single data files (i.e. N1C1W1_A, N1C1W1_B, ... N1C1W1_T were put
into the data file N1C1W1). The instances chosen are all from the set 'Scholl 1', a set of 720,
uniformly distributed instances with n between 50 and 500. The capacity c is between 100 and 150 (set ‘Scholl 1’).
The names of the data files are self-explanatory : Nx is the number of items (we have N1 = 50, N2 = 100, n3 = 120),
Cx is the capacity of the bins (C1 = 100, C2 = 120, C3 = 150), etc.
   
   
Please note that, according to Eilon and Christofides reported computational results, the 
branch-and-bound (BnB) algorithm implemented here can only solve small-size instances in a
reasonable time. Having that issue in mind, a timeout functionality was added to the BnB.
