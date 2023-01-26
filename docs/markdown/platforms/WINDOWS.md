Building GridPACK on Windows is not for the faint of heart. 

## System Preparation

The following are available as native Windows applications and can just be
installed in the normal Windows way:
* [Visual C++](https://www.visualstudio.com) in some form.
** The "free" ones, Visual Studio Express and Community, should work.
** To get the 64-bit compilers for Visual Studio 2015, you need to run the x64
   Native Tools Command Prompt (you may need to go into the Visual Studios 2015
   folder and then into the Visual Studio Tools subfolder and then into the Windows
   Desktop Command Prompts folder under that to find it). Once you have opened an
   x64 Native Command Prompt, cd to C:\Program Files (x86)\Microsoft Visual Studio
   14.0\VC and type `vcvarsall amd64`. More information about this can be found
   [here](https://msdn.microsoft.com/en-us/library/x4d2c09s.aspx).
* [Windows SDK](https://www.microsoft.com/en-us/download/details.aspx?id=8279)
  Usually, this is just installed with Visual Studio
* [CMake](https://cmake.org/download/) is required. During the installation
  process, the dialog will ask if you want CMake added to PATH variable. For
  some reason, the default is not to include CMake. You should change the
  selection to add CMake to the PATH for all users.
* A minimal [Cygwin](https://cygwin.com/) installation is necessary. Cygwin
  is required to build, and test for, [PETSc](https://www.mcs.anl.gov/petsc/). If
  the system has a Cygwin installation that is being used, it would probably be
  best to make another installation for GridPACK installation exclusively.  If
  `gcc` and `cmake` are already installed in Cygwin, you
  can rerun the installer and deselect these packages (they can be found under
  the developer category).
** Minimum required packages:
*** Base
*** Python
*** Make
** Rename `/usr/bin/link.exe` so it does not interfere with Windows `LINK.EXE` 
** If you want to use a Cygwin shell to build and/or debug GridPACK applications:
*** Do not install a compiler set.
*** Do not install CMake. 
* [https://msdn.microsoft.com/en-us/library/windows/desktop/bb524831%28v=vs.85%29.aspx
* Microsoft MPI: This appears to be the only modern implementation
  available for Windows.  In the past, the OpenMPI and MPICH implementations
  were available for Windows, but no more. You will need to download both the
  `.msi` and `.exe` files.
* Some software to unpack `.zip`, `.gz`, and `tar`
  archives.  Commands to unpack all of these archives are available with Cygwin.
  Windows can handle  `.zip` archives natively. 
* (optional) [MS-MPI Debugger * Extension](https://www.microsoft.com/en-us/download/details.aspx?id=48215)
  useful for debugging problem with parallel
  programs.

## Build Required Libraries 

In these instructions, everything is done from the Command Prompt.  Open a
64-bit Visual Studio command prompt, which should be available from the Start
Menu. If your version of Visual Studio has more than one type of Command
Prompt, use the developer prompt. You should bring up Command Prompts as
administrator, otherwise Visual Studio has trouble finding the common tools
folder. To do this, mouse over the Command Prompt field and right-click. A
menu should show up with "Run as Administrator" as one of the options. Type,
or copy and paste, the commands below into that command prompt window.
The caret `^` character is the line continuation character for Command
Prompt window.
Choose a place to install libraries. `C:\GridPACK` is used in this
case.  Avoid a path with spaces or special characters in it.  It's convenient
to set an environment variable to hold this path.  `GridPACKDir` is
used here.  
Make a folder in the install directory in which source archives can be
unpacked and the builds can be performed.  `C:\GridPACK\src` is used
here.  Again, avoid a path with spaces or special characters in it.
Prepare VS/CMake to use MPI. Enter the command
```
    set msmpi
```
The response should be something like 
```
    MSMPI_BIN=C:\Program Files\Microsoft MPI\Bin\
    MSMPI_INC=C:\Program Files (x86)\Microsoft SDKs\MPI\Include\
    MSMPI_LIB32=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x86\
    MSMPI_LIB64=C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\
```

## Boost

Download source from [here](https://sourceforge.net/projects/boost/files/) and
unpack in `%GridPACKDir%\src`
Remember to `set msmpi` as described above
Bootstrap the boost.build, e.g.,
```
    .\bootstrap.bat
```
Add `using mpi ;` to `project-config.jam`
As you might imagine, Boost is not built on Windows with MPI often.
Consequently, it's a little out of date. The file `mpi.jam` needs to be
modified. It's in a different place depending on the version. In Boost 1.61 it
is in `./tools/build/src/tools/mpi.jam`. Modify `mpi.jam` to
make it work using the following patch:
     247a248,250
     >     local win_ms_mpi_sdk = "C:\\Program Files (x86)\\Microsoft SDKs\\MPI"
     >     ;
     >     local win_ms_mpi = "C:\\Program Files\\Microsoft MPI" ;
     >
     249,251c252,254
     <     local cluster_pack_path_native = "C:\\Program Files\\Microsoft Compute Cluster Pack" ;
     <     local cluster_pack_path = [ path.make $(cluster_pack_path_native) ] ;
     <     if [ GLOB $(cluster_pack_path_native)\\Include : mpi.h ]
     ---
     >     # local cluster_pack_path_native = "C:\\Program Files\\Microsoft >     Compute Cluster Pack" ;
     >     # local cluster_pack_path = [ path.make $(cluster_pack_path_native) ]
     >     ;
     >     if [ GLOB $(win_ms_mpi_sdk)\\Include : mpi.h ]
     255c258
     <         ECHO "Found Microsoft Compute Cluster Pack: $(cluster_pack_path_native)" ;
     ---
     >         ECHO "Found Microsoft Compute Cluster Pack: $(win_ms_mpi_sdk)" ;
     260,262c263,265
     <       options = <include>$(cluster_pack_path)/Include
     <
<address-model>64:<library-path>$(cluster_pack_path)/Lib/amd64
     <                 <library-path>$(cluster_pack_path)/Lib/i386
     ---
     >       options = <include>$(win_ms_mpi_sdk)/Include
     >                 <address-model>64:<library-path>$(win_ms_mpi_sdk)/Lib/x64
     >                 <library-path>$(win_ms_mpi_sdk)/Lib/x86
     268c271
     <       .mpirun = "\"$(cluster_pack_path_native)\\Bin\\mpiexec.exe"\" ;
     ---
     >       .mpirun = "\"$(win_ms_mpi)\\Bin\\mpiexec.exe"\" ;

If you edit this file by hand, note that in addition to replacing the fields
`cluster_pack_path` and `cluster_pack_path_native`, you also need
to replace the values
`amd64` and `i386`.

* Configure, build, and install. The `python` library is not needed;
  `iostreams` is skipped unless some additional libraries are installed;
  `graph_parallel` fails to build (see
  [this ticket](https://svn.boost.org/trac/boost/ticket/11908)).  There are two
  options for building:  
* Option 1: just build what GridPACK requires (header-only libraries are still
installed)
```
  .\b2
    --prefix=%GridPACKDir%
    --with-mpi
    --with-serialization
    --with-random
    --with-filesystem
    --with-system
    --build-type=complete
    threading=single
    address-model=64
    link=static runtime-link=shared
    install
```
* Option 2: build everything (except python -- I'm not sure what's required for
that) 
```
  .\b2
    --prefix=%GridPACKDir%
    --without-python
    --build-type=complete
    threading=single
    address-model=64
    link=static runtime-link=shared
    install
```

Notes:

* 1.58.0 works with VS 2013
* 1.61.0 works with VS 2010, but requires update 5 for VS 2013 (see
  [this * ticket](https://svn.boost.org/trac/boost/ticket/11885)) 
* Boost decorates the library names with the compiler version, so explicitly
  specify the compiler and make sure that's the compiler you use for GridPACK.
  If multiple compilers are available, it is possible to force the boost
  configuration to pick a specific one by adding an option like
  `-toolchain=msvc-11.0` to `\.b2`.
  Some references for these instructions:
** A pretty complete set of instructions for
[build Boost on Windows with MPI](http://stackoverflow.com/questions/26147564/how-to-build-boost-mpi-for-ms-mpi-with-visual-studio-2012)
** [Some additional details](http://stackoverflow.com/questions/9433311/error-in-building-boost-mpi-in-msvc-2010/32635378#32635378).
** [Official Boost build instructions](http://www.boost.org/doc/libs/1_61_0/more/getting_started/windows.html)
** A pertinent Boost [ticket](https://svn.boost.org/trac/boost/ticket/11908).

## Algebra Libraries

Ordinarily, several packages used by GridPACK can be built by having PETSc
download and build them as part of the PETSc build. Unfortunately, while this
works quite well on the Linux OS, it does not work in many cases for Windows.
Consequently, these packages must be downloaded and built separately. Note that
the builds below all assume that you are calling the cmake command in a separate
directory, below the main library directory.

### BLAS/LAPACK (CLAPACK)

Some implementation of BLAS/LAPACK is required for PETSc and some other
libraries. The implementation described here was chosen because it does not
require a Fortran compiler.  It is apparently *really* slow.  It will probably
be necessary to install Intel compilers and MKL in order to get improvement in
speed.  The Windows port of CLAPACK is described
[here](http://icl.cs.utk.edu/lapack-for-windows/clapack/).
* Get the source * [here](http://icl.cs.utk.edu/lapack-for-windows/clapack/clapack-3.2.1-CMAKE.tgz).
  This is a gzip'd tar archive.  This can be unpacked by
  [7-zip](http://www.7-zip.org/download.html), if available. Since Cygwin is
  required to build and use PETSc, the command 
```
   tar vxzf clapack-3.2.1-CMAKE.tgz
```
  can be executed in a Cygwin terminal to unpack the archive in the current
  directory
* Configure, build, and install
```
    cd clapack-3.2.1-CMAKE
    mkdir build
    cd build
    cmake -Wdev
        -G "Visual Studio 10 2010 Win64"
        -D CMAKE_INSTALL_PREFIX:PATH="%GridPACKDir%"
        ..
    cmake --build . --config Release
    cmake --build . --target install --config Release
```
  Depending on which version of Visual Studio you are using, the argument for -G
  will change. You can get a listing of options by typing `cmake -G`.

### SuiteSparse

* For use with PETSc
* Get the source for the Windows port from
  [here](https://github.com/jlblancoc/suitesparse-metis-for-windows).
  Needs BLAS/LAPACK that's consistent with PETSc, so try to just use the one
  PETSc likes. Note that you can verify that both SuiteSparse and PETSc are
  using the same BLAS/LAPACK libraries by explicitly specifying them in both the
  SuiteSparse configuration (via the configuration variables
  `SUITESPARSE_CUSTOM_BLAS_LIB` and
  `SUITESPARSE_CUSTOM_LAPACK_LIB`) and the PETSc configuration (via the
  configuration option `--with-blas-lapack-lib`).
* Edit `CMakeLists.txt` in the root directory by commenting out the
  following lines to remove the library suffix:
```
  ## get POSTFIX for lib install dir
  #if(CMAKE_SIZEOF_VOID_P MATCHES "8")
  #  set(LIB_POSTFIX "64" CACHE STRING "suffix for 32/64 inst dir placement")
  #else()
    set(LIB_POSTFIX "" CACHE STRING "suffix for 32/64 inst dir placement")
  #endif()
  mark_as_advanced(LIB_POSTFIX)
```
* Configure, build, and install
```
  cmake -Wdev
      -G "Visual Studio 10 2010 Win64"
      -D BUILD_METIS:BOOL=NO
      -D SUITESPARSE_INSTALL_PREFIX:PATH="%GridPACKDir%"
      -D SUITESPARSE_USE_CUSTOM_BLAS_LAPACK_LIBS:BOOL=ON
      -D SUITESPARSE_CUSTOM_BLAS_LIB:PATH=%prefix%\lib\blas.lib
      -D SUITESPARSE_CUSTOM_LAPACK_LIB:PATH=%prefix%\lib\lapack.lib
      -D CMAKE_INSTALL_PREFIX:PATH="%GridPACKDir%"
      ..
  cmake --build . --config Release
  cmake --build . --target install --config Release
```


== ParMETIS ==

* Get the official distribution from
  [here](http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz).
  The source distribution needs to modified to 
  Have CMake look for MPI in the correct way
  Edit `metis.h` in the `parmetis-4.0.3/metis/include` directory
  so that `REALTYPEWIDTH` is 64.
  Make sure `metis.h` is installed
  Here is the patch for the `CMakeLists.txt` file located in the root
directory:
```
  diff -r -u parmetis-4.0.3/CMakeLists.txt parmetis-4.0.3.fixed/CMakeLists.txt
  --- parmetis-4.0.3/CMakeLists.txt     2013-03-30 09:24:50.000000000 -0700
  +++ parmetis-4.0.3.fixed/CMakeLists.txt       2016-06-30 11:32:38.691121400
-0700
  @@ -15,6 +15,7 @@
   #   message(FATAL_ERROR "mpi is not found")
   # endif()
   # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
  +find_package(MPI REQUIRED)

   # Prepare libraries.
   if(SHARED)
  @@ -33,6 +34,7 @@
   include_directories(${METIS_PATH}/include)

   # List of directories that cmake will look for CMakeLists.txt
  +add_subdirectory(${METIS_PATH}/include)
   add_subdirectory(${METIS_PATH}/libmetis ${CMAKE_BINARY_DIR}/libmetis)
   add_subdirectory(include)
   add_subdirectory(libparmetis)
  Only in parmetis-4.0.3.fixed/: CMakeLists.txt~
  diff -r -u parmetis-4.0.3/metis/include/metis.h
parmetis-4.0.3.fixed/metis/include/metis.h
  --- parmetis-4.0.3/metis/include/metis.h      2013-03-30 09:24:50.000000000
-0700
  +++ parmetis-4.0.3.fixed/metis/include/metis.h        2016-06-30
11:08:21.791160800 -0700
  @@ -40,7 +40,7 @@
      32 : single precission floating point (float)
      64 : double precission floating point (double)
   --------------------------------------------------------------------------*/
  -#define REALTYPEWIDTH 32
  +#define REALTYPEWIDTH 64
```

* Configure and build
```
  cmake
      -G "Visual Studio 10 2010 Win64"
      -D BUILD_SHARED_LIBS:BOOL=NO
      -D METIS_INSTALL:BOOL=YES
      -D CMAKE_INSTALL_PREFIX:PATH="C:\GridPACK"
      ..
  cmake --build . --config Release
  cmake --build . --target install --config Release
  copy ..\metis\include\metis.h %prefix%\include
  copy libmetis\Release\metis.lib %prefix%\lib
```

You still need to copy `metis.h` and `metis.lib` to the
installation directory. For later versions of Visual Studio (e.g. 2015), you may
get errors that look like

```
   C:\Program Files (x86)\Windows Kits\10\Include\10.0.10240.0\ucrt\math.h(520):
   error C2059: syntax error: '('
[C:\GridPACK\src\parmetis-4.0.3\build\libmetis\metis.vcxproj]
```

These can be fixed by going into the `parmetis-4.0.3/metis/GKlib`
directory and commenting out the macro defining `rint(x)` at the bottom
of `gk_arch.h`

```
   #ifdef __MSC__
   /* MSC does not have rint() function */
   /*
   #define rint(x) ((int)((x)+0.5))
   */
   /* MSC does not have INFINITY defined */
   #ifndef INFINITY
   #define INFINITY FLT_MAX
   #endif
   #endif
```

### PETSc
These instructions mostly follow
[the discussion here](https://github.com/INMOST-DEV/INMOST/wiki/0206-Compilation-PETSc-Windows)
as well as the instructions for building
[Suitesparse on Windows with CMake](https://github.com/jlblancoc/suitesparse-metis-for-windows)

PETSc must be built within a Cygwin shell.  In this case, Cygwin was installed
in `C:\cygwin64`
As mentioned above, PETSc refuses to fetch and build many external packages
when configured on Windows.  They need to be built individually.  See
instructions above for the following:
* BLAS/LAPACK
* ParMETIS
* SparseSuite
Start a VS Command Prompt
* Remember to set the Microsoft MPI environment
  set msmpi
* Within the VS Command Prompt, start a minimal Cygwin terminal and get it ready
  c:\cygwin64\bin\mintty.exe
* Within the Cygwin shell, review the contents of the `PATH` environment
  variable. Take out anything that is not directly related to the PETSc build.
* Then, make sure Cygwin commands are available: 
```
  export PATH="/usr/bin:$PATH"
```
* Configure and build for real values
```
  prefix="/cygdrive/c/GridPACK"
  ./configure \
    PETSC_ARCH=mswin-cxx-real-opt \
    --with-cc="win32fe cl" \
    --with-clanguage=c++ \
    --with-c++-support=1 \
    --download-f2cblaslapack=0  \
    --with-blas-lapack-lib=[${prefix}/lib/lapack.lib,${prefix}/lib/blas.lib,${prefix}/lib/libf2c.lib]
\
    --download-superlu_dist=0 \
    --download-metis=0 \
    --with-metis=1 \
    --with-metis-include=${prefix}/include \
    --with-metis-lib=[${prefix}/lib/metis.lib] \
    --download-parmetis=0 \
    --with-parmetis=1 \
    --with-parmetis-include=${prefix}/include \
    --with-parmetis-lib=[${prefix}/lib/parmetis.lib] \
    --download-suitesparse=0 \
    --with-suitesparse=1 \
    --with-suitesparse-include="${prefix}/include" \
    --with-suitesparse-lib=[${prefix}/lib/libumfpack.lib,${prefix}/lib/libamd.lib,${prefix}/lib/libbtf.lib,${prefix}/lib/libcamd.lib,${prefix}/lib/libccolamd.lib,${prefix}/lib/libcolamd.lib,${prefix}/lib/libcholmod.lib,${prefix}/lib/libcxsparse.lib,${prefix}/lib/libklu.lib,${prefix}/lib/libspqr.lib,${prefix}/lib/libldl.lib,${prefix}/lib/suitesparseconfig.lib]
\
    --with-c-support=0 \
    --with-fortran=0 \
    --with-fc=0 \
    --with-precision=double \
    --with-scalar-type=real \
    --with-mpi-include=/cygdrive/c/Program\ Files\ \(x86\)/Microsoft\
SDKs/MPI/Include \
    --with-mpi-lib=['/cygdrive/c/Program\ Files\ \(x86\)/Microsoft\
SDKs/MPI/Lib/x64/msmpi.lib'] \
    --with-mpi-mpiexec=/cygdrive/c/Program\ Files/Microsoft\ MPI/Bin/mpiexec.exe
\
    --with-debugging=0 \
    --with-windows-graphics=0 \
    --with-x11=0 \
    --CFLAGS='-O2 -MD -wd4996' \
    --CXXFLAGS='-O2 -MD -wd4996 -EHsc'
  make
  make test
```
The `-EHsc` flag is to enable C++ exceptions. Change
`--with-scalar-type=real` to `--with-scalar-type=complex` to build
PETSc using complex numbers (you may also want to change the `PETSC_ARCH`
variable from `mswin-cxx-real-opt` to `mswin-cxx-complex-opt`).
The complex build appears to have a few problems. If the build stops, keep
typing `make` until it completes.
* You may have problems using DOS-type directory names as the arguments to
 options such as `--with-mpi-include`, particularly if they contain
 spaces. These can be converted to Linux-type directory names by using the
 cygpath command inside the Cygwin window. Type `cygpath -d "C:\Program
* Files (x86)\Microsoft SDKs\MPI\Include"` to get
 `C:\PROGRA~2\MIA713~1\MPI\Include` and use this name for any configure
 options that have a path listing.

We have had difficulties building PETSc using Visual Studio 2015. For PETSc
3.6, we encountered too many compiler errors to get it to build but with PETSc
3.7 we were able to successfully build by modifying the file
PETSC/src/mat/order/amd/amd.c. Change

```
  tval = (PetscBool)Control[AMD_AGGRESSIVE];
```

to
```
  /* tval = (PetscBool)Control[AMD_AGGRESSIVE]; */
  if (Control[AMD_AGGRESSIVE]) tval = PETSC_TRUE;
  else tval = PETSC_FALSE;
```

* Verify installation by building the tests in GridPACK
`sandbox/petsc-cmake`:
```
  set path=%path%;c:\cygwin64\bin
  set CFLAGS="/D _ITERATOR_DEBUG_LEVEL=0"
  set CXXFLAGS="/D _ITERATOR_DEBUG_LEVEL=0"
  cmake -Wdev
    -D BOOST_ROOT:PATH=C:\GridPACK
    -D Boost_USE_STATIC_LIBS:BOOL=ON
    -D BOOST_INCLUDEDIR=C:\GridPACK\include\boost-1_61
    -D PETSC_DIR:PATH="C:\GridPACK\src\petsc-3.6.4"
    -D PETSC_ARCH:STRING='mswin-cxx-complex-opt'
    -G "Visual Studio 10 2010 Win64"
    ..
  cmake --build . --config Release
```

### Global Arrays

* Get the source from the trunk Subversion repository (see the
  [home page](http://www.emsl.pnl.gov/docs/global home page)):
* Note that you need to use the CMake build for GA.
* Make a directory in which to build and change to that directory
```
   mkdir ga-svn\build
   cd ga-svn\build
```
* Configure, build, and install
```
  set CFLAGS="/D _ITERATOR_DEBUG_LEVEL=0"
  set CXXFLAGS="/D _ITERATOR_DEBUG_LEVEL=0"
  cmake -Wdev --debug-trycompile
    -G "Visual Studio 10 2010 Win64"
    -D ENABLE_BLAS:BOOL=No
    -D ENABLE_FORTRAN:BOOL=No
    -D ENABLE_CXX:BOOL=Yes
    -D GA_RUNTIME:STRING=MPI_TS
    -D CMAKE_INSTALL_PREFIX:PATH="%GridPACKDir%"
    ..
  cmake --build . --config Release
  cmake --build . --config Release --target install
```

* Check using Global Arrays tests

* Check w/ GridPACK sandbox. In the `sandbox/ga` directory, add a build
  directory and cd into it. Configure the GA test with
```
  cmake -Wdev --debug-trycompile
    -G "Visual Studio 10 2010 Win64"
    -D BOOST_ROOT:PATH="%GridPACKDir%"
    -D Boost_USE_STATIC_LIBS:BOOL=ON
    -D BOOST_INCLUDEDIR="%GridPACKDir%\include\boost-1_61"
    -D GA_DIR:PATH="%GridPACKDir"
    ..
```

= GridPACK =

* Check out the [Windows GridPACK fork](https://github.com/wperkins/GridPACK/tree/windoze)
```
  git clone -b windoze https://github.com/wperkins/GridPACK.git gridpack
  cd gridpack
  git submodule init
  git submodule update
```
* Make a build directory and change into it
```
  mkdir build
  cd build
```
* Configure
```
    set path=%path%;c:\cygwin64\bin
    set prefix="C:\GridPACK"
    set CFLAGS="/D _ITERATOR_DEBUG_LEVEL=0"
    set CXXFLAGS="/D _ITERATOR_DEBUG_LEVEL=0"
    cmake -Wdev --debug-trycompile
      -G "Visual Studio 10 2010 Win64"
      -D USE_PROGRESS_RANKS:BOOL=OFF
      -D BOOST_ROOT:PATH=C:\GridPACK
      -D Boost_USE_STATIC_LIBS:BOOL=ON
      -D Boost_USE_DEBUG_RUNTIME:BOOL=OFF
      -D BOOST_INCLUDEDIR=C:\GridPACK\include\boost-1_61
      -D PETSC_DIR:PATH="C:\GridPACK\src\petsc-3.6.4"
      -D PETSC_ARCH:STRING='mswin-cxx-complex-opt'
      -D GA_DIR:PATH='C:\GridPACK\ga-svn'
      -D GA_TEST_RUNS:BOOL=YES
      -D PARMETIS_DIR:PATH=C:\GridPACK
      -D MPIEXEC_MAX_NUMPROCS:STRING="2"
      -D GRIDPACK_TEST_TIMEOUT:STRING=60
      ..
```
  Note: this is the contents of `gridpack/example_configuration.bat`
* Build
```
  cmake --build . --config Release
```
* Install
```
  cmake --build . --config Release --target install
```
