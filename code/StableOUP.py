import scipy.stats as stats
import numpy as np
import pandas as pd
from multiprocessing import Pool
from ast import literal_eval
import sys

class StableOUP:

	def gamma_equilibrium(alpha):
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

		samplesize = 1
		for k in K:
			for s in range(0, samplesize):

				sample_k = []

				for n in range(0,5):

					dL = stats.levy_stable.rvs(alpha = alpha, beta = 0, loc = 0, scale = 1, random_state = (5+s)*n, size = N)
					dtdL = dL*(dt)**(1/alpha)

					X = np.zeros(N)
					Xtemp = Xzero

					for j in range(1, N):
						Linc = 0.1*dtdL[j]
						Xtemp = Xtemp - dt*k*Xtemp + Linc
						X[j-1] = Xtemp


					sample_k.append(X[-70:-1])



				paras = stats.levy_stable.fit(np.hstack(sample_k))
				gamma_sample[i] = paras[3]
				run[i] = s
				K_in[i] = k

				i += 1

				df = pd.DataFrame({'gamma_sample': gamma_sample,
						   'k': K_in, 'run': run}).to_csv("data/new/gamma_lin_eq_alpha" + str(alpha) + "_rev.csv")
						   
	def gamma_neq(alpha):
		print(alpha)
		s = 1
		K = np.round(np.arange(5,0, -0.0001), 3)
		k = 5
		window = 300 #every so many timesteps are we calculating gamma
		
		Xzero = 0.5

		dt = 0.005
		N = len(K)
		L = 10000 #length spinup
		
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
						    'k': K_in, 'run': run, 'timestep': timestep, 'sample': S}).to_csv("data/new/gamma_lin_neq_alpha" + str(alpha) + "_" + str(s) + ".csv") #write a file for each
						    
	def trajectories_neq(alpha):
		# problem parameters
	
		K = np.round(np.arange(100,0, -0.001), 3)
		
		k = 5	
		Xzero = 0.5

		T = 5
		N = 2**10
		dt = 1/N

		R = 4
		Dt = R*dt
		L = int(N/R)*T
		t = np.arange(dt, T+dt, dt)

		#allocate space to solutions
		l = 500

		gamma_sample = np.zeros(l)
		run = np.zeros(l)
		K_in = np.zeros(l)

		i = 0

		samplesize = 1
		
		df = pd.DataFrame({'template': np.zeros(L + len(K)), 'K': np.concatenate((np.repeat(k, repeats = L), np.array(K)))}).to_csv("data/new/trajectories_lin_neq_" + str(alpha) + ".csv", index=False)

		for s in range(0, samplesize):
			print(s)

			
			dL = stats.levy_stable.rvs(alpha = alpha, beta = 0, loc = 0, scale = 1, random_state = i, size = N*T + len(K)*R)
			dtdL = dL*(dt)**(1/alpha)
		
			X = np.zeros(L + len(K))
			print(len(X))
			print(L + len(K))
			Xtemp = Xzero
			for j in range(1, L):
				Linc = 0.1*sum(dtdL[R*(j-1):R*j])
				Xtemp = Xtemp - k*Xtemp*Dt + Linc
				X[j-1] = Xtemp
				
			for j in range(0, len(K)):	
				j += L
				Linc = 0.1*sum(dtdL[R*(j-1):R*j])
				Xtemp = Xtemp - K[j-L]*Xtemp*Dt + Linc
				X[j-1] = Xtemp
					
					
			traj = pd.read_csv("data/new/trajectories_lin_neq_" + str(alpha) + ".csv")
			traj['i' + str(s) +  '_k' + str(k)] = X
			traj.to_csv("data/new/trajectories_lin_neq_" + str(alpha) + ".csv", index=False)
			


			i += 1
			


    


