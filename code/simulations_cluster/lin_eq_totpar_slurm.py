import scipy.stats as stats
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from multiprocessing import Pool
import os
from ast import literal_eval
import sys

# problem parameters


def simulate(alpha, k):
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
	for s in range(0, samplesize):

		sample_k = []

		for n in range(0,5):

			dL = stats.levy_stable.rvs(alpha = alpha, beta = 0, loc = 0, scale = 1, random_state = (5+s)*n, size = N)
			dtdL = dL*(dt)**(1/alpha)

			X = np.zeros(N)
			Xtemp = Xzero

			for j in range(1, N):
				Linc = 0.1*dtdL[j]
				Xtemp = Xtemp - (dt*k*Xtemp)/(1 + dt*abs(k*Xtemp)) + Linc
				X[j-1] = Xtemp


			sample_k.append(X[-70:-1])



		paras = stats.levy_stable.fit(np.hstack(sample_k))
		gamma_sample[i] = paras[3]
		run[i] = s
		K_in[i] = k

		i += 1

		df = pd.DataFrame({'gamma_sample': gamma_sample, 'k': K_in, 'run': run}).to_csv("gamma_lin_eq_alpha" + str(alpha) + "_k" + str(k) + ".csv") #write a file for each K


var1=literal_eval(sys.argv[1])
var2=literal_eval(sys.argv[2])


simulate(var1, var2)


