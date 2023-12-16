import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import scipy.stats as stats
import time
import datetime
from multiprocessing import Pool
import os
from ast import literal_eval
import sys



def simulate(alpha, s):  #define a function that calculates one trajectory of one alpha	
	print(alpha)
	print(s)
	K = np.round(np.arange(5,0, -0.0001), 3)
	k = 5
	window = 300 #every so many timesteps are we calculating gamma, large for noe to test
	
	Xzero = 0.5

	dt = 0.005
	N = len(K)
	L = 10000 #length spinup -> now same length as equilibrium simulation. it does not cost much and only then can we be sure we are in eq
	
	t = [k]*L

	t.extend(K)
		

	#allocate space to solutions
	l = int(round(N/window, 0)) + 100
	gamma_sample = np.zeros(l)
	run = np.zeros(l)
	timestep = np.zeros(l)
	K_in = np.zeros(l)
	S = np.zeros(l)
	S = s

	i = 0 #start iteration
	
	dL = stats.levy_stable.rvs(alpha = alpha, beta = 0, loc = 0, scale = 1, random_state = s, size = N + L)
	dtdL = dL*(dt)**(1/alpha)   

	X = np.zeros(N + L)
	Xtemp = Xzero

	for j in range(1, L): #spin simulation up until equilibrium
		Linc = 0.1*dtdL[j]
		Xtemp = Xtemp - dt*k*Xtemp + Linc
		X[j-1] = Xtemp
						
	for j in range(0, N):  #start estimation gamma every window timestep
		Linc = 0.1*dtdL[j]
		Xtemp = Xtemp - dt*K[j]*Xtemp + Linc
		X[L + j-1] = Xtemp

		if j%window == 0:
			paras = stats.levy_stable.fit(X[L + j-300:L + j - 1])
			gamma_sample[i] = paras[3]
			run[i] = s
			timestep[i] = j
			K_in[i] = K[j]
			i += 1

			df = pd.DataFrame({'gamma_sample': gamma_sample,
					    'k': K_in, 'run': run, 'timestep': timestep, 'sample': S}).to_csv("gamma_lin_neq_alpha" + str(alpha) + "_" + str(s) + ".csv") #write a file for each trajectory
					    
					    
			
var1=literal_eval(sys.argv[1])
var2=literal_eval(sys.argv[2])

simulate(var1, var2)
				   
