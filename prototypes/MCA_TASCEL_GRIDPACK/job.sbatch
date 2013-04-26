#!/bin/csh
#SBATCH -A MCA 
#SBATCH -t 00:29:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J runca
#SBATCH --exclusive
#SBATCH -p short
##SBATCH --partition=gpu
##SBATCH --partition=pal
##SBATCH --partition=xl
#SBATCH -o runca.out.%j 
#SBATCH -e runca.err.%j
#SBATCH --mail-user=yousu.chen@pnl.gov

# First make sure the module commands are available.
#source /etc/profile.d/modules.csh 

#Then, set up your environment you wish to run in with module commands.
module purge
module load precision/i4
module load intel 
module load openmpi 
module load mkl/11.1
#module load pathscale

#module load openmpi

#Next unlimit system resources, and set any other environment variables you need.
unlimit

#Now it is extremely useful to record the modules you have loaded, your limit settings, 
#your current environment variables and the dynamically load libraries that your executable 
#is linked against in your job output file.
echo
echo "loaded modules"
echo
module list >& _modules.lis_
cat _modules.lis_
/bin/rm -f _modules.lis_
echo
echo limits
echo
limit
echo
echo "Environment Variables"
echo
printenv
echo
echo "ldd output"
echo
ldd runca

#Now you put your parallel launch command in.
#For each different parallel executable you launch we recommend
#adding a corresponding ldd command to verify that the environment
#that is loaded corresponds to the environment the executable was built in.

#srun --mpi=none runca temp 
#mpirun --mca btl openib,self,sm --mca btl_openib_receive_queues X,128,256,192,128:X,2048,256,128,32:X,12288,256,128,32:X,65536,256,128,32 runca 
mpirun runca

