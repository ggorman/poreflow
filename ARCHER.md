Tutorial for running on ARCHER (www.archer.ac.uk)
==========================

If you have not already done so then read the [basic tutorial]{https://github.com/ggorman/poreflow/blob/master/tutorial.md} first.

~/.bashrc file setup
--------------------
I recommend appending these lines to your ~/.bashrc file if they are not already there:
```bash
export MODULEPATH=$PWD/modules:/work/e319/shared/modules:$MODULEPATH
export WORK=/work/n02/n02/$USER/ # This path will vary depending on what project you are on. The important part is that you are operating from /work and not /home.

export PACKAGE_VERSION="dev"
export INSTANT_CACHE_DIR=$WORK/.instant
```

Download and compile
--------------------

```bash
mkdir -p $WORK/projects
cd $WORK/projects
git clone https://github.com/ggorman/poreflow.git
cd poreflow
```
As it turns out we do not have to build any code on ARCHER - we are only going to use the python code that uses Dolfin.

Prepare data sample
-------------------
This should already have been completed following the instructions from the [basic tutorial]{https://github.com/ggorman/poreflow/blob/master/tutorial.md}. So create a suitable directory on ARCHER (within /work of course) and scp across the dolfin .xml data file.

Running a simulation
--------------------
When using ARCHER (or indeed any cluster/supercomputer) you typically do not run interactively like on a workstation. Instead you use a batch queuing system - such as PBS in this case. Here is an example PBS file for submitting a job to the queue.
```bash
cat > mysimulation.pbs << EOF
#
# Parallel script produced by bolt
#        Resource: ARCHER (Cray XC30 (24-core per node))
#    Batch system: PBSPro_select
#
# bolt is written by EPCC (http://www.epcc.ed.ac.uk)
#
#PBS -l select=4
#PBS -N pore
# The resources account this simulation is charged to.
#PBS -A n02-ic1
# The maximum wall time this simulation will run to.
#PBS -l walltime=8:00:00

# Switch to the directory where you submitted the pbs script to the queue.
cd $PBS_O_WORKDIR

# Load the fenics environment lovingly maintained on ARCHER by Chris Richardson and Garth Wells (thanks guys).
module load fenics/dev-64bit

# Use an alternative gcc version to avoid some pesky bugs.
module switch gcc gcc/4.9.0

# -n specifies the total number of processing elements (pes - i.e. cores)
# -N number of pes that we are using per node
# -S number of pes we are using per NUMA region.
# -d number of cpu's per pe
# In general you just want to specify -l select=x to specify the number of compute nodes you want, and then update the -n option on the next line to 24*number_of_nodes. You can leave the other options as they are.
aprun -n 96 -N 24 -S 12 -d 1 python $WORK/projects/poreflow/stokes-dolfin-snes.py Berea.xml
EOF
To submit this "job" to the queue execute the command:
```bash
qsub mysimulation.pbs
```
You monitor whether it is sitting in the queue or actually running using the command:
```bash
qstat -u $USER
```
If you do not get anything back that means that your job has completed. At this point you should look at your output files:
```bash
ls -ltr
```
You will see files of the pattern pore.o* and pore.e* which are the standard output and standard error of your job respectively. 

Possible problems
-----------------
Sometimes the "instant" run compiler fails if it runs in parallel. The way around this problem is to first run a small problem in serial so that all the runtime compiling completes. Afterward you can run in parallel without issue.

