import pandas as pd
import os
import glob


for a in [2, 1.8, 1.5, 1.3]:
	files = glob.glob('*benchmark_gammaX_alpha' + str(a) + '*.{}'.format('csv'))

	print(files)


	pd.concat([pd.read_csv(f) for f in files ], ignore_index=True).to_csv("benchmark_gammaX_alpha" + str(a) + ".csv")
	
