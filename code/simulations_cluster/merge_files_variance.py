import pandas as pd
import os
import glob


for a in [2, 1.5]:
	files = glob.glob('*convergance_alpha' + str(a) + '*.{}'.format('csv'))

	print(files)


	pd.concat([pd.read_csv(f) for f in files ], ignore_index=True).to_csv("convergance_alpha" + str(a) + ".csv")
	
