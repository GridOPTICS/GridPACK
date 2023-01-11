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

Example scripts for configuring and building Boost can be found on the links
below.

* [Mac Yosemite ](../DUMMY.md)
* [Redhat Linux Workstation](../DUMMY.md)
* [CentOS 6](../DUMMY.md)
* [Redhat Linux Cluster](../RC_CLUSTER.md)
