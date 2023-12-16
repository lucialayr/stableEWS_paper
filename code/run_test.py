from StableFold import *
from StableOUP import *


if __name__ == '__main__':
    alpha_values = [2, 1, 1.5, 0.5]
    
    for a in alpha_values:

	    StableOUP.gamma_equilibrium(alpha=a)
	    #StableOUP.gamma_neq(alpha=a)
	    StableOUP.trajectories_neq(alpha=a)
	   

	    StableFold.gamma_equilibrium(alpha=a)
	    #StableFold.gamma_neq(alpha=a)
	    StableFold.trajectories_neq(alpha=a)
