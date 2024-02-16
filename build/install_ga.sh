#! /bin/bash

ga_version="${GA_VERSION:-5.8}"

echo "--- Installing Global Arrays ${ga_version} ---"

# download
echo "Downloading Global Arrays"
wget \
  "https://github.com/GlobalArrays/ga/releases/download/v${ga_version}/ga-${ga_version}.tar.gz" \
  -O ga.tar.gz \
  --quiet

# unpack
echo "Unpacking Global Arrays"
tar -xf ga.tar.gz && rm -f ga.tar.gz

# remove version from dir name
mv "ga-${ga_version}" ga

pushd ga || exit

# build
echo "Configuring Global Arrays"
./configure \
  --with-mpi-ts \
  --disable-f77 \
  --without-blas \
  --enable-cxx \
  --enable-i4 \
  --prefix="${PWD}/install_for_gridpack" \
  --enable-shared

# install
echo "Installing Global Arrays"
make -j "${MAKE_JOBS:-$(nproc)}" install

popd || exit

echo "Global Arrays installation complete"
