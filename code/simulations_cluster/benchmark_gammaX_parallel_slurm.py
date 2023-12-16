import scipy.stats as stats
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from ast import literal_eval
import sys



def simulate(a, nt, w):
	k = 1
	
	Xzero = 0.5

	dt = 0.005
	N = 10**3
	T = N/dt
	t = np.arange(dt, T+dt, dt)
	
	samplesize = 100

	#allocate space to solutions
	l = samplesize + 100

	gamma_sample = np.zeros(l)
	run = np.zeros(l)
	W_in = np.zeros(l)
	T_in = np.zeros(l)
	alpha_in = np.zeros(l)

	i = 0

	

	for s in range(0, samplesize):
		sample_k = []		

		for n in range(1, nt + 1):
			dL = stats.levy_stable.rvs(alpha = a, beta = 0, loc = 0, scale = 1, random_state = 5 + s*n, size = N)
			dtdL = dL*(dt)**(1/a)

			X = np.zeros(N)
			Xtemp = Xzero

			for j in range(1, N):
				Linc = 0.1*dtdL[j]
				Xtemp = Xtemp - dt*k*Xtemp + Linc
				X[j-1] = Xtemp

			sample_k.append(X[-w:-1])


		paras = stats.levy_stable.fit(np.hstack(sample_k))
		gamma_sample[i] = paras[3]
		run[i] = s
		W_in[i] = w
		T_in[i] = nt
		alpha_in[i] = a


		df = pd.DataFrame({'gamma_sample': gamma_sample, 'alpha': alpha_in, 'run': run,
				    'W_in': W_in, 'T_in': T_in}).to_csv("benchmark_gammaX_alpha" + str(a) + "_w" + str(w) + "s" +  str(nt) + ".csv")
				    
		i += 1
					   
var1=literal_eval(sys.argv[1])
var2=literal_eval(sys.argv[2])
var3=literal_eval(sys.argv[3])

simulate(var1, var2, var3)
