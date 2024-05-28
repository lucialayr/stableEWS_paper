#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --error=%j.err
#SBATCH --output=%j.log
#SBATCH --job-name=gamma_eq
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


for a in 2 1.8 1.5 1.3

	do
	
		for k in 100 10 1 0.1 0.001
		
		do
		
		srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python lin_eq_totpar_slurm.py ${a} ${k} &
		srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python nol_eq_totpar_slurm.py ${a} ${k} &


		done

	done

wait

srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python merge_files_eq.py &

wait


