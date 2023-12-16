import scipy.stats as stats
import time
import datetime
from multiprocessing import Pool
import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import math

from multiprocessing import Pool
import os
from ast import literal_eval
import sys


def simulate(alpha):
	K = [100, 10, 1, 0.1, 0.01]
	Xzero = 0.5


	dt = 0.005
	N = 10**4
	T = N/dt
	
	t = np.arange(dt, T+dt, dt)

	#allocate space to solutions
	l = 500


	gamma_sample = np.zeros(l)
	run = np.zeros(l)
	K_in = np.zeros(l)

	i = 0

	samplesize = 100
	for k in K:
		
		X_critical = - np.sqrt(k) - k*0.1 #threshhold after which we consider trajectory to be escaped beyond unstable FP
		
		for s in range(0, samplesize):

			sample_k = []

			for n in range(0,5):

				dL = stats.levy_stable.rvs(alpha = alpha, beta = 0, loc = 0, scale = 1, random_state = (5+i)*n, size = N)
				dtdL = dL*(dt)**(1/alpha)

				X = np.zeros(N)
				Xtemp = Xzero

				for j in range(1, N):
					Linc = 0.1*dtdL[j]
					Xtemp = Xtemp + dt*(k - Xtemp**2) + Linc
					X[j-1] = Xtemp
					if Xtemp < X_critical:
						print("overflow encountered at timestep " + str(j))
						break

				sample_k.append(X[j-75:j-5])
				
			paras = stats.levy_stable.fit(np.hstack(sample_k))
			gamma_sample[i] = paras[3]
			run[i] = s
			K_in[i] = k

					

			df = pd.DataFrame({'gamma_sample': gamma_sample,
					   'k': K_in, 'run': run}).to_csv("gamma_nol_eq_alpha" + str(alpha) + "_rev.csv")
			i += 1
			
var1=literal_eval(sys.argv[1])

number_nodes = len(var1)

if __name__ == '__main__':
    with Pool(number_nodes) as p:
        print(p.map(simulate, var1))
				   
