#!/bin/bash
# compiling the program
mpifort -O2 cholesky_mpi_block_test.f90 -o cholesky_mpi_block_test
# Entering the maximum number of processes, the amount of experiments and the block width
echo "Enter the maximum number of processes"
read max_nproc
echo "Enter the amount of experiments"
read nexp
echo "Enter the block width"
read nb
# Running the program in the following loop
for (( exp=1; exp<=$nexp; exp++ ))
do
   echo "Experiment: $exp"
   for (( nproc=1; nproc<=$max_nproc; nproc++ ))
   do
      echo "Number of processes: $nproc"
      filename=../data/$exp.$nproc.txt
      echo "# n; Execution time, sec" > $filename
      ns=(100 1100 2500 5000 10000)
      for n in ${ns[@]}
      do
         echo "n: $n"
         # keeping time data
         echo "$n $nb" | mpirun -np $nproc ./cholesky_mpi_block_test >> $filename
      done
   done
done
