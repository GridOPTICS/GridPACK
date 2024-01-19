#! /bin/bash

# installs necessary packages for a given distribution

set -xeuo pipefail

distribution=$(
  # shellcheck source=/dev/null
  source /etc/os-release
  echo "$ID"
)

case $distribution in
debian | ubuntu)
  apt update &&
    apt upgrade -y &&
    apt install -y \
      wget build-essential git python3.11 python3-pip libopenmpi-dev cmake pkg-config
  ;;
fedora | rhel | centos | rocky)
  dnf update -y &&
    dnf install -y \
      wget @development git python3.11 python3-pip openmpi-devel cmake pkgconf
  ;;
*)
  echo "$distribution not supported"
  exit 1
  ;;
esac