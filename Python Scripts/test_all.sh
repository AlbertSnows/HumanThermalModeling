#!/bin/bash -l
##$ -l h_rt=0:02:00

for i in 1 2
  do
     sbatch --constraint=elves --time=180 --mem-per-cpu=10240 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$i --nodes=1 run_unmodified.sh
done
sbatch --constraint=elves --time=180 --mem-per-cpu=6144 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=4 --nodes=1 run_unmodified.sh


for j in 1 2
  do
     sbatch --constraint=elves --time=180 --mem-per-cpu=10240 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$j --nodes=1 run_serial_optimized.sh
done
sbatch --constraint=elves --time=180 --mem-per-cpu=6144 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=4 --nodes=1 run_serial_optimized.sh


for k in 1 2
  do
     sbatch --constraint=elves --time=180 --mem-per-cpu=10240 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$k --nodes=1 run_multithread.sh
done
sbatch --constraint=elves --time=180 --mem-per-cpu=6144 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=4 --nodes=1 run_multithread.sh


for l in 1 2
  do
     sbatch --constraint=elves --time=180 --mem-per-cpu=10240 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=$l --nodes=1 run_multiprocess.sh
done
sbatch --constraint=elves --time=180 --mem-per-cpu=6144 --partition=ksu-gen-reserved.q,batch.q --ntasks-per-node=4 --nodes=1 run_multiprocess.sh