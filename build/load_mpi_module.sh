#! /bin/bash

# source this script to load the mpi module on RHEL distributions

distribution=$(
  source /etc/os-release
  echo "$ID"
)

case $distribution in
  fedora | rhel | centos | rocky)
    echo "Loading mpi module"
    source /etc/profile.d/modules.sh
    module load "mpi/openmpi-$(arch)"
    ;;
esac
