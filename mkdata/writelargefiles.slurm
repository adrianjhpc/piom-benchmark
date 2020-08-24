#!/bin/bash --login
#SBATCH -D /home/fh01/fh01/adrianj/piom-benchmark/mkdata
#SBATCH -o /home/fh01/fh01/adrianj/piom-benchmark/mkdata/mkdata.%A.out
#SBATCH -e /home/fh01/fh01/adrianj/piom-benchmark/mkdata/mkdata.%A.err
#SBATCH -J mkdata
#SBATCH --nodes=8
#SBATCH --ntasks=512
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:0

XPROCS=32
YPROCS=16
BARRIER=0

EXECUTABLELOCATION=/home/fh01/fh01/adrianj/piom-benchmark/mkdata

cd /lustre/hdd/fh01/fh01/adrianj/piom-benchmark/data/striped

echo "max stripes"
#srun $EXECUTABLELOCATION/parallelcmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
#srun $EXECUTABLELOCATION/parallelfmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
srun $EXECUTABLELOCATION/parallelcmkinfileindivid 640000 144000 $XPROCS $YPROCS $BARRIER

cd /lustre/hdd/fh01/fh01/adrianj/piom-benchmark/data/unstriped

echo "1 stripe"
#srun $EXECUTABLELOCATION/parallelcmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
#srun $EXECUTABLELOCATION/parallelfmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
srun $EXECUTABLELOCATION/parallelcmkinfileindivid 640000 144000 $XPROCS $YPROCS $BARRIER

