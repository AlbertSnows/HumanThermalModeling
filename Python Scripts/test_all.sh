#!/bin/bash -l
##$ -l h_rt=0:02:00

for i in 1 2 4 8 16
  do
     sbatch --constraint=dwarves --time=360 --mem-per-cpu=1024 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$i --nodes=1 run_unmodified.sh
done

for j in 1 2 4 8 16
  do
     sbatch --constraint=dwarves --time=360 --mem-per-cpu=1024 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$j --nodes=1 run_serial_optimized.sh
done

for k in 1 2 4 8 16
  do
     sbatch --constraint=dwarves --time=360 --mem-per-cpu=1024 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$k --nodes=1 run_multithread.sh
done

for l in 1 2 4 8 16
  do
     sbatch --constraint=dwarves --time=360 --mem-per-cpu=1024 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$l --nodes=1 run_multiprocess.sh
done
