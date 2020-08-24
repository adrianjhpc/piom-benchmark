#!/bin/bash --login

#!/bin/bash
#SBATCH -D /home/nx01/nx01/adrianj/piom-benchmark/mkdata
#SBATCH -o /home/nx01/nx01/adrianj/piom-benchmark/mkdata/mkdata.%A.out
#SBATCH -e /home/nx01/nx01/adrianj/piom-benchmark/mkdata/mkdata.%A.err
#SBATCH -J mkdata
#SBATCH --nodes=8
#SBATCH --ntasks=384
#SBATCH --cpus-per-task=1
#SBATCH --nvram-options=1LM:1000
#SBATCH --time=2:00:0

export PSM2_MULTI_EP=1
export PSM2_MULTIRAIL=1
export PSM2_MULTIRAIL_MAP=0:1,1:1

export FI_PROVIDER=psm2
export FI_PSM2_LAZY_CONN=0
export PSM2_DEVICES=self,hfi,shm


NPROCS=384
XPROCS=32
YPROCS=12
BARRIER=0

EXECUTABLELOCATION=/home/nx01/nx01/adrianj/piom-benchmark/mkdata

cd /home/nx01/nx01/adrianj/piom-benchmark/data/striped

echo "max stripes"
#srun $EXECUTABLELOCATION/parallelcmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
#srun $EXECUTABLELOCATION/parallelfmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
srun $EXECUTABLELOCATION/parallelcmkinfileindivid 640000 144000 $XPROCS $YPROCS $BARRIER

cd /home/nx01/nx01/adrianj/piom-benchmark/data/unstriped

echo "1 stripe"
#srun $EXECUTABLELOCATION/parallelcmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
#srun $EXECUTABLELOCATION/parallelfmkinfile 640000 144000 $XPROCS $YPROCS $BARRIER
srun $EXECUTABLELOCATION/parallelcmkinfileindivid 640000 144000 $XPROCS $YPROCS $BARRIER

