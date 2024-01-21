#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --error=%j.err
#SBATCH --output=%j.log
#SBATCH --job-name=benchmark_gamma
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<emailadress>
#SBATCH --time=100:00:00
#SBATCH --clusters=htls
#SBATCH --partition=htls_batch
#SBATCH --reservation=htls_users
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --get-user-env
#SBATCH --export=NONE

module load slurm_setup
module load anaconda3

source activate mypython38

srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python benchmark_gammaL.py &

for a in 2 1.8 1.5 1.3
do
	for nt in 1 #3 5 10
	do
	
		for w in 10 50 70 100
		do
		
		 srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python benchmark_gammaX_parallel_slurm.py ${a} ${nt} ${w} &
		 
		done
		
	done
	
done

wait

srun -c 1 --mpi=none --mem-per-cpu=1800MB --exclusive -n 1 -N 1 python merge_files_neq.py &
wait



