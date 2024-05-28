import pandas as pd
import os
import glob


for a in [2, 1.8, 1.5, 1.3]:
	linear_files = glob.glob('*lin_eq_alpha' + str(a) + '*.{}'.format('csv'))

	print(linear_files)


	pd.concat([pd.read_csv(f) for f in linear_files ], ignore_index=True).to_csv("gamma_lin_eq_alpha" + str(a) + ".csv")
	
	nol_files = glob.glob('*nol_eq_alpha' + str(a) + '*.{}'.format('csv'))

	print(nol_files)


	pd.concat([pd.read_csv(f) for f in nol_files ], ignore_index=True).to_csv("gamma_nol_eq_alpha" + str(a) + ".csv")
