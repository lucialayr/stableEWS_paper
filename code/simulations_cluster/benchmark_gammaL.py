import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import scipy.stats as stats
import datetime


#benchmark parameter estimation of free alpha-stable process
alpha = [2, 1.8, 1.5, 1.3]

samplesize = 200

T = 5
N = 2**10
dt = 1/N
#allocate space for output

fitted_gamma = np.zeros(samplesize*len(alpha)+100)
fitted_gamma_scaled = np.zeros(len(fitted_gamma))
true_gamma = np.zeros(len(fitted_gamma))
true_alpha = np.zeros(len(fitted_gamma))
simulation = np.zeros(len(fitted_gamma))
i = 0

data_length = 5*70 #what we use in paper as well

for a in alpha:
	for s in range(0,samplesize):
		dL = stats.levy_stable.rvs(alpha=a, beta=0, scale = 1, loc = 0, size=5*70)
		dt_dL = dL*dt**(1/a)
		
		paras_dL = stats.levy_stable.fit(dL)
		paras_dt_dL = stats.levy_stable.fit(dt_dL)
		
		simulation[i] = s
		true_alpha[i] = a
		true_gamma[i] = 1
		fitted_gamma[i] = paras_dL[3]
		fitted_gamma_scaled[i] = paras_dt_dL[3]
		
		df = pd.DataFrame({'simulation': simulation, 'alpha':true_alpha, 'true_gamma': true_gamma, 'gamma': fitted_gamma, 'gamma_scaled': fitted_gamma_scaled,}).to_csv("benchmark_gammaL.csv")
		i += 1
		
		print("first round was sucessful!")
				
				



