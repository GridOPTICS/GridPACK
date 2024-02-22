#! /bin/bash

# library to install packages required to build GridPACK on debian or RHEL based distros
# meant to be sourced

# load the mpi module on RHEL distributions
function load_mpi_module {
  local base_distro
  base_distro=$(get_base_distro)

  if ! test "$base_distro" = "rhel"; then
    return
  fi

  echo "Loading mpi module"

  # clear unbound variable check temporarily
  option_state=$([[ $- == *u* ]] && echo "-u" || echo "+u")
  set +u

  # source the script which contains the module command
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh

  # reset the unbound variable check setting
  set "$option_state"

  # load the mpi module
  module load "mpi/openmpi-$(arch)"
}

# detect if this is a "debian" or "rhel" based linux distribution
# usage:
#   if test "$(get_base_distro)" = "debian"; then
#      apt-get update
#      apt-get upgrade --yes
#   fi
function get_base_distro {
  distribution=$(
    # shellcheck disable=SC1091
    source /etc/os-release
    echo "$ID"
  )

  case $distribution in
  debian | ubuntu)
    echo "debian"
    ;;
  fedora | rhel | centos | rocky)
    echo "rhel"
    ;;
  *)
    echo "unknown"
    ;;
  esac
}

# install packages for RHEL
function install_rhel_packages {
  dnf upgrade --assumeyes --verbose

  # use epel and crb repos
  # https://wiki.rockylinux.org/rocky/repo/#notes-on-epel
  dnf install epel-release --assumeyes
  crb enable

  # install required packages
  dnf install --assumeyes \
    wget @development git python3.11 python3-pip openmpi-devel cmake pkgconf \
    python3-mpi4py-openmpi
}

# install packages for debian
function install_debian_packages {
  apt-get update
  apt-get upgrade --yes

  # install required packages
  apt-get install --yes \
    wget build-essential git python3.11 python3-pip libopenmpi-dev cmake pkg-config python3-mpi4py
}

function install_packages {
  local base_distro
  base_distro=$(get_base_distro)

  if test "$base_distro" = "debian"; then
    install_debian_packages
  elif test "$base_distro" = "rhel"; then
    install_rhel_packages
  fi
}
