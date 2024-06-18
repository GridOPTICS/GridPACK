# Quick guide for GridPACK-EMT
## Installation
GridPACK depends on three external packages (Boost, Global Arrays (GA), and PETSc). This section guides through the installation of these packages via a script (`install_gridpack_deps.sh`) we've developed.

### Prerequisites
Before beginning with the GridPACK installation process, we need to have some prerequisite libraries installed. These libraries are the fundamental building blocks of parallel applications and their build process.  

- Make sure you have a functional MPI installation. You can use either OpenMPI or MPICH. MPI is a prerequisite for GridPACK and its dependencies.
- Install CMake ver. 3.18.1 or higher. This is required for PETSc dependencies. You could also have PETSc install CMake. This can be done by adding `--download-cmake` to the PETSc configure. (line 118 in install_gridpack_deps.sh).
- You also need to have functional `Blas/Lapack` libraries that can be found by the installer. If you want PETSc to install Blas/lapack then add `--download-f2cblaslapack=1` to the PETSc configure line (line 118 in install_gridpack_deps.sh).

### Install GridPACK dependencies
1. Go to the root GridPACK directory.
```
cd $GRIDPACK_ROOT_DIR
```
Here, `$GRIDPACK_ROOT_DIR` is the root GridPACK directory


3. Assuming you have installed MPI and CMake, source `install_gridpack_deps.sh` to install all GridPACK dependencies. This will install all the GridPACK dependencies Boost (ver. 1.81), Global Arrays (ver. 5.8), and PETSc (ver. 3.20.1) in `$GRIDPACK_ROOT_DIR/external-dependencies`.
```
source install_gridpack_deps.sh
```

If you wish to install the libraries at another location then set `GP_EXT_DEPS` environment variable to the location you want to install the libraries. For example,
```
export GP_EXT_DEPS=/usr/local
```
will install all the external packages in `usr/local`

4. Once the dependencies get installed, it is best to set some environmental variables and update library paths so that the linker is always able to find the libraries (this is only to do if you have installed the dependencies in a non-standard location)
```
# Root GridPACK directory                                                                                        
#export GRIDPACK_ROOT_DIR=<location_of_GridPACK_directory>                      
# GridPACK external dependencies                                                                            
#export GP_EXT_DEPS=${GRIDPACK_ROOT_DIR}/external-dependencies                                                     
# GridPACK install directory                                                                                       
#export GRIDPACK_DIR=${GRIDPACK_ROOT_DIR}/src/install                                                              
# update LD_LIBRARY_PATH so that boost,ga, and petsc are on it                                                     
#export LD_LIBRARY_PATH=${GP_EXT_DEPS}/boost_1_78_0/install_for_gridpack/lib:${GP_EXT_DEPS}/ga-5.8/install_for_gri\
dpack/lib:${GP_EXT_DEPS}/petsc/install_for_gridpack/lib:${LD_LIBRARY_PATH}
```
### Installation errors
Here are some of the installation errors we've seen on different systems.
1. <b>PETSc cannot find Blas/Lapack:</b> This is because PETSc cannot find a functional Blas/Lapack. The error looks like this...
```
--------------------------------------------------------------------------------------------
  Could not find a functional BLAS. Run with --with-blas-lib=<lib> to indicate the library
  containing BLAS.
  Or --download-f2cblaslapack=1 to have one automatically downloaded and installed
*********************************************************************************************

Building PETSc 3.20.1
makefile:24: /home/abhy245/software/GridPACK/external-dependencies/petsc/build-dir/lib/petsc/conf/petscrules: No such file or directory
make[1]: *** No rule to make target '/home/abhy245/software/GridPACK/external-dependencies/petsc/build-dir/lib/petsc/conf/petscrules'.  Stop.
...
...
```
To fix this, add `--download-f2cblaslapack=1` to the PETSc configure line (line 118 in install_gridpack_deps.sh) and re-run the installation script.

2. <b>PETSc cannot find CUDA libraries</b>

If your system has NVIDIA CUDA libraries then PETSc may throw the following warning when completing its installation. This is a harmless warning. <b> You can ignore it. </b>

```
...
Running PETSc check examples to verify correct installation
Using PETSC_DIR=/people/abhy245/software/GridPACK/external-dependencies/petsc and PETSC_ARCH=build-dir
Possible error running C/C++ src/snes/tutorials/ex19 with 1 MPI process
See https://petsc.org/release/faq/
--------------------------------------------------------------------------
The library attempted to open the following supporting CUDA libraries,
but each of them failed.  CUDA-aware support is disabled.
libcuda.so.1: cannot open shared object file: No such file or directory
libcuda.dylib: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.so.1: cannot open shared object file: No such file or directory
/usr/lib64/libcuda.dylib: cannot open shared object file: No such file or directory
If you are not interested in CUDA-aware support, then run with
--mca mpi_cuda_support 0 to suppress this message.  If you are interested
in CUDA-aware support, then try setting LD_LIBRARY_PATH to the location
of libcuda.so.1 to get passed this issue.
--------------------------------------------------------------------------
lid velocity = 0.0016, prandtl # = 1., grashof # = 1.
Number of SNES iterations = 2
Possible error running C/C++ src/snes/tutorials/ex19 with 2 MPI processes
See https://petsc.org/release/faq/
--------------------------------------------------------------------------
```

<b>Note:</b> If the error is only with PETSc install then you can have the installation script just install PETSc. This can be done by setting the `install_boost` flag (line 7) and `install_ga` flag (line 8) to `false`. This way, only PETSc installation will be executed.

## Installing GridPACK
Once the dependencies are installed, you can install the GridPACK library by running
```
source install_gridpack.sh
```
This will build and install the GridPACK library. The build files will be in `$GRIDPACK_ROOT_DIR/src/build` and the install files will be in `$GRIDPACK_ROOT_DIR/src/install`

## Executing GridPACK-EMT
1. Go to the `bin` folder in GridPACK installation directory 
` cd $GRIDPACK_ROOT_DIR/src/install/bin`. This is the directory from where you'll run the EMT simulation.
2. Each GridPACK-EMT simulation run needs four input files.
- A PETS configuration file where all PETSc configuration options are specified. This file needs to have the name `.petscrc` for PETSc to recognize it.
- The input network file in PTI raw format.
- The input dynamic data file with the extension `.dyr`
- A GridPACK configuration xml file that specifies the configuration for GridPACK and the EMT simulation parameters

GridPACK provides some sample data sets and the `.petscrc` that we can use. So, lets copy these files over to the current directory in the next step.
3. Next, we'll copy over some input files needed for executing the EMT application. The `$GRIDPACK_ROOT_DIR/src/applications/data_sets/emt` directory has sample EMT data sets for a 2-bus system, 9-bus system, Kundur two-area network, and 39-bus system. Let's copy over the input data files from the `two-bus` directory.
```
cp $GRIDPACK_ROOT_DIR/src/applications/data_sets/emt/two-bus/input_2bus.xml .
cp $GRIDPACK_ROOT_DIR/src/applications/data_sets/emt/two-bus/case2.dyr .
cp $GRIDPACK_ROOT_DIR/src/applications/data_sets/emt/two-bus/case2mod.raw .
```
Additionally, copy-over the PETSc configuration file
```
cp $GRIDPACK_ROOT_DIR/src/applications/data_sets/petscoptions/emt/.petscrc .
```
Your current directory should now have the files for the two-bus system and the `.petscrc` file. Note, you would not be able to see the `.petscrc` file if you simply do `ls`. You will need to do `ls -a`.
### GridPACK configuration file
The `input_2bus.xml` is the GridPACK configuration that controls the input for the EMT simulation. Let's take a look at the main elements of this xml file. Note, the file has some other elements for controlling other features but we will only focus on those that are needed for EMT.
1. Setting the raw file: The network raw file is set with the tag `<networkConfiguration>`. In the `input_2bus.xml` file the network file is set to
`<networkConfiguration> case2mod.raw </networkConfiguration>`. Recall `case2mod.raw` is the file we copied. GridPACK supports various version of PTI raw data files including v23, v33, v34, v35. `<networkConfiguration>` tag assumes the file format is v23. If you have files in other versions then add an `_vXX` to the tag where `XX` is the version number. For e.g., to read a v34 file use `<networkConfiguration_v34>` tag.
2. Setting the dyr file: The dyr file is set with the tag `<generatorParameters>`. In the `input_2bus.xml` file the dyr file is set to 
`<generatorParameters> case2.dyr </generatorParameters>`
3. Setting the simulation time and timestep: The simulation time and time-step is set with the tags `<simulationTime>` and `timeStep`, respectively.
4. Scenarios/Disturbances/Events: GridPACK currently supports two different types of disturbances. These events need to be within the `<Events>` tag. 
-  Bus faults: Bus faults can be single-phase or three-phase. A three-phase temporary fault can be set for the simulation with the following tags.
```
<BusFault>
  <begin>0.5</begin>
  <end>0.51</end>
  <bus>5</bus>
  <type>ThreePhase</type>
  <Ron>1e-3</Ron>
  <Rgnd>1e-2</Rgnd>
</BusFault>
```
Here, `<begin>` is the fault begin time, `<end>` is the fault end time, `<bus>` is the bus number where fault is incident, `<type>` should be `ThreePhase`, `Ron` and `Rgnd` are the per unit fault and ground resistances, respectively.

Similarly, a single-phase to ground fault is specified with
```
<BusFault>
  <begin>0.5</begin>
  <end>0.51</end>
  <bus>5</bus>
  <type>SLG</type>
  <phase>A</phase>
  <Ron>1e-3</Ron>
  <Rgnd>1e-2</Rgnd>
</BusFault>
```
Here, the `<type>` should be `SLG` and the `<phase>` where fault is applied should be specified.

- A generator trip example configuration is shown below

```
<GenTrip>
  <bus>5</bus>
  <id>1</id>
  <phase>A</phase>
  <time>1e-3</time>
</BusFault>
```
Here, `<bus>` and `<id>` specify the generator bus number and id, respectively, and `<time>` is the time at which the generator is tripped.

5. Output monitoring: The output from EMT simulation can be saved to output csv file. GridPACK currently supports saving a subset of generator output (generator terminal voltage, generator real power, angle,and per unit speed)
```
  <Monitors>
    <Generator>
      <bus> 1 </bus>
      <id> 1 </id>
    </Generator>
  </Monitors>
```
and three-phase instantaneous bus voltages
```
  <Monitors>
    <Bus>
      <bus> 1 </bus>
    </Bus>
  </Monitors>
```
 By default, the output gets saved in `emtoutput.csv`. The output can be set with the following tags. You can multiple `<Generator>` and `<Bus>` monitors.

 ### Example simulation run
 The EMT application is run with the following command.
 ```
 ./emt.x <gridpack_configuration_file>
 ``` 
 where `<gridpack_configuration_file>` is the xml configuration file.
 
 For the two-bus system data files we copied, running the emt simulation will result in the following output.
 ```
 ./emt.x input_2bus.xml
 ```
 At the end of the simulation run, you should see a file name `emtoutput.csv` where the output for Generator at bus 1 is saved.





 

