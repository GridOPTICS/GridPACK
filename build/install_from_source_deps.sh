#! /bin/bash

# installs GridPACK dependencies to the current directory

# bash options:
# - xtrace: print each command before executing it
# - errexit: exit on error
# - nounset: treat unset variables as errors
# - pipefail: treat whole pipeline as errored if any commands within error
# https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
set -o xtrace -o errexit -o nounset -o pipefail

# get the parent directory of this script
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

echo "Installing GridPACK dependencies"
date

bash "${script_dir}/install_boost.sh"
bash "${script_dir}/install_ga.sh"
bash "${script_dir}/install_petsc.sh"

echo "GridPACK dependency installation complete"
