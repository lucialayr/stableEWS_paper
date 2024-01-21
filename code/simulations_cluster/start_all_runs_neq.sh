#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --error=%j.err
#SBATCH --output=%j.log
#SBATCH --job-name=simulate_more_alphas_parallel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<enteremail>
#SBATCH --time=100:00:00
#SBATCH --clusters=htls
#SBATCH --partition=htls_batch
#SBATCH --reservation=htls_users
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=3
#SBATCH --get-user-env
#SBATCH --export=NONE

module load slurm_setup
module load anaconda3

source activate mypython38


for a in 2 1.8 1.5 1.3

	do
	
		for i in {1..10}
		
		do
		
		srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python lin_neq_totpar_slurm.py ${a} ${i} &
		srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python nol_neq_totpar_slurm.py ${a} ${i} &


		done

	done

wait

srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python merge_files_neq.py &

wait


