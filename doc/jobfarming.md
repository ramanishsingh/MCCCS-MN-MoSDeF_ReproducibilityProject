# Job farming in MCCCS-MN
MCCCS-MN utilizes a job-farming module to perform large numbers of MC simulation in parallel in a single job. Each simulation in the farmed job is run independently so that the job farming module can be readily scaled to tens of thousands of CPUs.

## Compilation
When compiling the code, use the flag `-DJOB_FARMING=ON` or turn on the `JOB_FARMING` switch in CMake GUI to enable job farming. Each CPU core performs a single simulation so we uuse `-DUSE_MPI=ON` and `-DUSE_OPENMP=OFF`.

### Compiling on MSI Mangi (AMD Epyc CPUs):
Below is a working configuration for job farming on MSI for AMD CPUs:
1. Load the GCC compiler suite, OpenMPI library, and cmake. While we can compile MCCCS-MN for individual simulations using Intel compilers, **Job farming runs only work on GCC and OpenMPI.**
```bash
module load gcc/8.2.0 ompi/4.0.0/gnu-8.2.0-centos7 cmake
```
2. Generate the makefile for MCCCS-MN using cmake. You can also use CMake GUI with `ccmake /path/to/MCCCS-MN`.
```bash
export FC=gfortran && export CC=gcc && export CXX=g++
cmake -DJOB_FARMING=ON -DUSE_MPI=ON -DUSE_OPENMP=OFF /path/to/MCCCS-MN
make -j16
```
## Job script and Directory structure

Below is a typical directory structure of running a farmed job:
```
/
|--job.sh
|
|--simulations
|  |--MFI-0-P0-T0
|  |  |...
|  |
|  |--MFI-0-P1-T0
|  |  |... 
|  |
|  |--MFI-0-P2-T0
|  |  |...
|  |   
|  |...
|  
|--jobs
   |--run1.txt
   |--(done.txt)
   |--(log.txt)
```
* The `simulations` directory and its subfolders contain the input files (`topmon.inp`, `fort.4`, etc.) for each simulation. The output files will also be generated in the same directory as input files.

* The `jobs/run1.txt` lists all the directories where each of the simulations is run. In other words, each worker core in the farmed job will `cd` into the directory listed in this file and runs a simulation just as the regular MCCCS-MN executable. The `run1.txt` file will look like (with 2 independent simulations per state point):
```
simulations/MFI-0-P0-T0/1
simulations/MFI-0-P0-T0/2
simulations/MFI-0-P1-T0/1
simulations/MFI-0-P1-T0/2
simulations/MFI-0-P2-T0/1
simulations/MFI-0-P2-T0/2
...
```
* After a farmed job has finished, two additional output files will be generated: `jobs/done.txt` and `jobs/log.txt`. The `done.txt` file records the simulations that has successfully finished. This file **must be removed before running the next batch of farmed job**. The `log.txt` records which core runs which simulation, for example, "`Rank 1234: MFI-0-P0-T0/1`". This information is useful when finding the simulation that causes error exit of the farmed job.

* In the job script, the MCCCS-MN executable is invoked with three arguments `1 1 1`:

`mpirun -np $N /path/to/topmon 1 1 1`

where `$N` is the number of simulations (number of lines in the `jobs/run1.txt` file) **plus one** for an additional master rank to dispatch the simulations.

## Input file settings
Simulations with different system compositions run at different speeds. Therefore in a farmed job we typically set a simulation **to run for a fixed time duration** instead of a fixed number of Monte Carlo cycles. This is done by adding an additional ``time_limit = $X`` config in the first section of `fort.4` file where `$X` is the simulation time limit in seconds. Also make sure that the number of cycles is large enough (e.g., 200000) to make sure the simulation uses the entire time limit.

Our current job farming module is implemented in the most straightforward way that the whole farmed job will exit even if a single simulation fails. Therefore simulations in a farmed job should be configured in a way **that the chance of the simulation error exits is minimized.** Below are the special input file settings that has been used for job farming simulations:
* In production runs, set `allow_cutoff_failure = 1`. This prevents the whole job from crashing when several of the simulation with the largest densities have simulation boxes smaller than the cutoff radius.
* Use small maximum displacement updates in early stages of equilibration, such as `iratio = 50` and `iratv = 50`. This is because the boxlength may change dramatically in early equilibration, and the simulation will error exit when the maximum translation displacement or maximum volume displacement is too large.

