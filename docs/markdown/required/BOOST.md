## Boost

The [Boost C++ Library](http://www.boost.org/) is used heavily throughout the
GridPACK framework, and a relatively recent version is required.  The GridPACK
configuration requires version 1.49 or later.  The
[Boost](http://www.boost.org/) installation must include
[Boost::MPI](http://www.boost.org/doc/libs/1_53_0/doc/html/mpi.html)
which must have been built with the same MPI compiler used for GridPACK.
Be aware that many Boost modules do not include MPI, so you may have to build
Boost on your own even if it is available as a module on your system.

To configure GridPACK to recognize Boost, one need only specify where
[Boost](http://www.boost.org/) is installed, like this

```
    -D BOOST_ROOT:STRING='/path/to/boost' \
```

Boost is tied quite closely to the latest features in C++ and problems can be
encountered if the version of Boost that you are using was released much later
than the compiler. Reverting to an earlier Boost version can sometimes eliminate
problems if you are having difficulties building it. The same is true for Boost
and CMake. If the CMake version was released earlier than the Boost version,
CMake may have problems identifying the libraries in Boost that it needs for
GridPACK. Again, going to an earlier version of Boost may fix these issues.

The Boost build can be tricky. Some clusters have Boost modules that can
potentially be used instead of building Boost on your own, but many modules are
not built with Boost::MPI. You will still need to specify the location of the
Boost directory, which can be found by using the command

```
    module show boost
```

This should tell you the location of `BOOST_ROOT`, which you can then use
in your GridPACK configuration script. If the GridPACK configuration does not
report that MPI was found, then you will need to get your system administrator
to rebuild Boost with the MPI libraries or build boost on your own. A successful
Boost configuration in GridPACK should report the results

```
    -- Checking Boost ...
    -- Boost version: 1.61.0
    -- Found the following Boost libraries:  
    --   mpi
    --   serialization
    --   random
    --   filesystem
    --   system
```

If you need to build Boost yourself, refer to the documentation on building
GridPACK on individual platforms for additional details on build Boost. If an
attempt to configure and build Boost fails, it usually is a good idea to fix the
build script and then remove the existing Boost directory and create a new one
by untarring the Boost tarball. Attempts to resume a failed Boost build after
fixing the build script are usually unsuccessful.

Download the Boost tarfile (e.g. boost_1_78_0.tar.gz). Older versions of Boost
can be found [here](https://www.boost.org/users/history/). The following script
should work with many versions of Boost. Information on downloading tar files
and creating scripts can be found [here](LINUX_BASICS.md#linux-basics).

A basic script for build Boost using GNU compilers and creating static libraries
is

```
    echo "using mpi ;" > ~/user-config.jam
    sh ./bootstrap.sh \
        --prefix="$PREFIX" \
        --without-icu \
        --with-toolset=gcc \
        --without-libraries=python,log
    ./b2 -a -d+2 link=static stage
    ./b2 -a -d+2 link=static install
    rm ~/user-config.jam
```

Run this script from the top level Boost directory. This should configure, build
and install the Boost libraries.

(In Boost 1.54,
[Boost.log](http://www.boost.org/doc/libs/1_54_0/libs/log/doc/html/index.html)
was added.  This uses some compiler capabilities not supported by the compilers
that come with older versions of RHEL/CentOS, so Boost.Log is
disabled.  Boost seems to work fine this way with older versions of RHEL.)

To build using the Intel compilers, substitute
`--with-toolset=intel-linux` for `--with-toolset=gcc`. You
may also run into problems with the name of the MPI wrapper for the C++
compiler. If it looks like configure is not finding `mpic++` then
replace the first line in the above script with

    echo "using mpi : /absolute/path/to/mpi/C++/wrapper ;" > ~/user-config.jam

Make sure you include the spaces around ":" and before ";".

Boost has a tendency to use cutting-edge features of the C++ compiler so it is a
good idea to use a compiler version that was released at the same time as the
Boost version you are working with. If you are having problems, you may have
better luck moving to an earlier version of Boost. If the Boost build fails, you
should delete the entire boost directory and start from scratch after making
corrections to your build script. Restarting a failed Boost build does not
appear to work in most instances.

If you want to use Intel compilers modify the `--with-toolset=gcc` line to
`--with-toolset=intel-linux`. For shared library builds, modify the two link lines
to

```
    ./b2 -a -d+2 link=shared stage
    ./b2 -a -d+2 link=shared install
```
