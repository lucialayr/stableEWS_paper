import scipy.stats as stats
import time
import datetime
from multiprocessing import Pool
import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from ast import literal_eval
import sys

# problem parameters

alpha = [2, 1.5]

def simulate(alpha, s):
	k = 5
	
	Xzero = 0.5

	dt = 0.005
	N = 5*10**4
	T = N/dt

	#allocate space to solutions
	steps = int(N/200 + 2)
	samplesize = 100
	l = steps + 100


	gamma_sample = np.zeros(l)
	var_sample = np.zeros(l)
	run = np.zeros(l)
	K_in = np.zeros(l)
	time_step = np.zeros(l)
	
	i = 0
	
	window = 300
	
	
	print("run " + str(s) + " for alpha " + str(alpha))

	dL = stats.levy_stable.rvs(alpha = alpha, beta = 0, loc = 0, scale = 1, random_state = (5+s), size = N)
	dtdL = dL*(dt)**(1/alpha)

	X = np.zeros(N)
	Xtemp = Xzero
	
	step_counter = 0

	for j in range(1, N):
		Linc = 0.1*dtdL[j]
		Xtemp = Xtemp - (dt*k*Xtemp)/(1 + dt*abs(k*Xtemp)) + Linc
		X[j-1] = Xtemp

		
		
		if j%window == 0:
			time_step[i] = j
			run[i] = s
			K_in[i] = k
			var_sample[i] = np.var(X[j-199:j])
			paras = stats.levy_stable.fit(X[j-199:j])
			gamma_sample[i] = paras[3]
	

			i += 1

			df = pd.DataFrame({'gamma_sample': gamma_sample, 'var': var_sample, 'time_step': time_step,
			   	'k': K_in, 'run': run}).to_csv("convergance_alpha" + str(alpha) + "_s" + str(s) + ".csv")
				
var1=literal_eval(sys.argv[1])
var2=literal_eval(sys.argv[2])

simulate(var1, var2)									   
