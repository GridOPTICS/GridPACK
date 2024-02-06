#! /bin/bash

# source this script to load the mpi module on RHEL distributions

# bash options:
# - xtrace: print each command before executing it
# - errexit: exit on error
# - nounset: treat unset variables as errors
# - pipefail: treat whole pipeline as errored if any commands within error
# https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
set -o xtrace -o errexit -o nounset -o pipefail

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
