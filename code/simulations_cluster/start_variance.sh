#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --error=%j.err
#SBATCH --output=%j.log
#SBATCH --job-name=variance_convergance
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lucia.layritz@tum.de
#SBATCH --time=100:00:00
#SBATCH --clusters=htls
#SBATCH --partition=htls_cm4
#SBATCH --reservation=htls_users
#SBATCH --ntasks=160
#SBATCH --ntasks-per-core=2
#SBATCH --get-user-env
#SBATCH --export=NONE


module load slurm_setup
module load anaconda3

source activate mypython38

for a in 2 1.5
do
	for i in {1..112}
	do
	
		
		
	 srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python variance_convergance_parallel.py ${a} ${i} &
		 

		
	done
	
wait
	
done



srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python merge_files_variance.py &

wait



