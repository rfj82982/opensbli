#!/bin/sh
#! Example SLURM job script for Iridis GPU
#! Last updated: Mon Aug 22 11:17:10 BST 2022
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#! SBATCH -J opensbli

#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<=nodes*12)
#SBATCH --ntasks=1
#! Specify the number of GPUs per node (between 1 and 4; must be 4 if nodes>1).
#! Note that the job submission script will enforce no more than 3 cores per GPU:
#! SBATCH --gres=gpu:1
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#! How much wallclock time will be required?
#SBATCH --time=00:10:00
#! What types of email messages do you wish to receive?
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#!SBATCH --no-requeue

#! Do not change:
#! SBATCH -p ampere

#! sbatch directives end here (put any additional directives above this line)

#! Notes:
#! Charging is determined by GPU number*walltime.

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
#. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
# module load rhel8/default-amp

#! Insert additional module load commands after this line if needed:
export OPS_COMPILER=intel
export NV_ARCH=Tesla

module load gcc/10.3.0
module load cuda/11.7
module load hdf5/1.12.2/intel/parallel

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/local/software/hdf5/1.12.2/gcc/parallel/lib:/local/software/cuda/11.7/lib:/local/software/gcc/10.3.0/lib
export HDF5_PATH=/local/software/hdf5/1.12.2/gcc/parallel/
#export I_MPI_DEBUG=1
#! Full path to application executable:
cd $SLURM_SUBMIT_DIR
# make opensbli_cuda

application="./opensbli_cuda"

#! Run options for the application:
#options="369 768 1024 OPS_TILING OPS_TILING_MAXDEPTH=6"
#options="OPS_FORCE_DECOMP_Y=4 OPS_FORCE_DECOMP_Z=2"
options=""
#options="0 0 0"
#options="369 768 1024"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 12:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
CMD="$application $options &> out.txt.$SLURM_JOB_ID"

#! Choose this for a MPI code using OpenMPI:
#CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options --mca mtl ^mxm --mca pml ^yalla --mca btl self,openib --mca btl_openib_want_cuda_gdr 1"

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

#! echo -e "JobID: $JOBID\n======"
#! echo "Time: `date`"
#! echo "Running on master node: `hostname`"
#! echo "Current directory: `pwd`"

#! if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        #! export NODEFILE=`generate_pbs_nodefile`
        #! cat $NODEFILE | uniq > machine.file.$JOBID
        #! echo -e "\nNodes allocated:\n================"
        #! echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
#! fi

#! echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

#! echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD





