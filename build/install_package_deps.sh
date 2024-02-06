#! /bin/bash

# installs necessary packages for a given distribution

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

  debian | ubuntu)

    apt-get update
    apt-get upgrade --yes

    # install required packages
    apt-get install --yes \
      wget build-essential git python3.11 python3-pip libopenmpi-dev cmake pkg-config python3-mpi4py
    ;;

  fedora | rhel | centos | rocky)

    dnf upgrade --assumeyes --verbose

    # use epel and crb repos
    # https://wiki.rockylinux.org/rocky/repo/#notes-on-epel
    dnf install epel-release --assumeyes
    crb enable

    # install required packages
    dnf install --assumeyes \
      wget @development git python3.11 python3-pip openmpi-devel cmake pkgconf \
      python3-mpi4py-openmpi
    ;;

  *)

    echo "$distribution not supported"
    exit 1
    ;;

esac
