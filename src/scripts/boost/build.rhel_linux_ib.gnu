#
# This script has been used to build Boost 1.65 on a Linux cluster using
# OpenMPI 1.8.3 and GNU 4.9.2 compilers. The cluster is using an Infiniband
# interconnect. Users should modify the argument of the --prefix option to
# point to the top level Boost directory where this script should be run.
#
echo "using mpi ;" > ~/user-config.jam
sh ./bootstrap.sh \
    --prefix="/pic/projects/gridpack/software_new/boost_1_65_0" \
    --without-icu \
    --with-toolset=gcc \
    --without-libraries=python
./b2 -a -d+2 link=static stage
./b2 -a -d+2 link=static install
rm ~/user-config.jam
