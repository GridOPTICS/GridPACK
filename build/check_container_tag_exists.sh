#! /bin/bash

# creates a file named container-tag-not-found if a tag does not exist in the container registry
# needs env vars CI_JOB_TOKEN, CI_API_V4_URL, CI_PROJECT_ID, IMAGE_TAG, CI_PROJECT_PATH_SLUG
# needs curl and jq


# bash options:
# - xtrace: print each command before executing it
# - errexit: exit on error
# - nounset: treat unset variables as errors
# - pipefail: treat whole pipeline as errored if any commands within error
# https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
set -o xtrace -o errexit -o nounset -o pipefail


# check that a command is available
check_installed() {
  if ! command -v "$1" &>/dev/null; then
    echo "please install '$1'". >&2
    exit 1
  fi
}


# submit a GET request to some endpoint on the project api
function get_container_registry_repos {
  # get endpoint or use empty string if not provided
  # https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html
  local endpoint=${1:-}

  # curl options:
  # - no-progress-meter: don't show progress bar but show errors
  # - location: follow redirects
  # - fail: return 22 on server errors
  # - header: send private token for authentication
  # https://docs.gitlab.com/ee/api/container_registry.html
  curl \
    --no-progress-meter \
    --fail \
    --location \
    --header "PRIVATE-TOKEN: ${CI_JOB_TOKEN:?}" \
    "${CI_API_V4_URL:?}/projects/${CI_PROJECT_ID:?}/registry/repositories${endpoint}"
}


check_installed jq
check_installed curl

# query the registry repo id
repo_id=$(
  # jq query:
  # - raw-output: output raw strings instead of json
  # - exit-status: return 1 if no results
  # - for each elements of the array
  # - where the registry repo path is the one we care about
  # - select the id prop
  # https://jqlang.github.io/jq/manual/
  get_container_registry_repos | jq --raw-output --exit-status ".[] | select(.path == \"${CI_PROJECT_PATH_SLUG}\") | .id"
)

# using that repo id check for a particular tag
# if command errors, create file named container-tag-not-found
get_container_registry_repos "/${repo_id}/tags/${IMAGE_TAG:?}" >/dev/null ||
  touch container-tag-not-found
