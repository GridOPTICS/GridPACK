#! /bin/bash

boost_version="${BOOST_VERSION:-1.78.0}"

echo "--- Installing Boost ${boost_version} ---"

# remove existing
rm -rf boost*

# download
echo "Downloading Boost"
wget \
  "https://boostorg.jfrog.io/artifactory/main/release/${boost_version}/source/boost_${boost_version//./_}.tar.gz" \
  -O boost.tar.gz \
  --quiet

# unpack
echo "Unpacking Boost"
tar -xf boost.tar.gz && rm -f boost.tar.gz

# remove version from dir name
mv "boost_${boost_version//./_}" boost

pushd boost || exit

# bootstrap
echo "Bootstrapping Boost"
./bootstrap.sh \
  --prefix=install_for_gridpack \
  --with-libraries=mpi,serialization,random,filesystem,system
echo 'using mpi ;' >>project-config.jam

# build
echo "Building Boost"
./b2 -a -d+2 link=shared stage

# install
echo "Installing Boost"
./b2 -a -d+2 link=shared install

popd || exit

echo "Boost installation complete"
