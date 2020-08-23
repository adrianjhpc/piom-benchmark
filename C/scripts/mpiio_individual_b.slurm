#!/bin/bash --login

#SBATCH -D /home/nx01/nx01/adrianj/piom-benchmark/C
#SBATCH -o /home/nx01/nx01/adrianj/piom-benchmark/C/indiv_b.%A.out
#SBATCH -e /home/nx01/nx01/adrianj/piom-benchmark/C/indiv_b.%A.err
#SBATCH -J indiv_b
#SBATCH --nodes=8
#SBATCH --tasks-per-node=48
#SBATCH --nvram-options=1LM:1000
#SBATCH --time=0:30:0

export PSM2_MULTI_EP=1
export PSM2_MULTIRAIL=1
export PSM2_MULTIRAIL_MAP=0:1,1:1

export FI_PROVIDER=psm2
export FI_PSM2_LAZY_CONN=0
export PSM2_DEVICES=self,hfi,shm

XPROCS=32
YPROCS=12
BARRIER=1

echo "1 stripe"
srun ./mpiio_individual /home/nx01/nx01/adrianj/piom-benchmark/data/unstriped/cinput640000x144000.dat 640000 144000 $XPROCS $YPROCS $BARRIER
echo "max stripes"
srun ./mpiio_individual /home/nx01/nx01/adrianj/piom-benchmark/data/striped/cinput640000x144000.dat 640000 144000 $XPROCS $YPROCS $BARRIER
