#! /bin/bash

# library of functions related to building the container image

# print a set of args as a comma-space seperated list
# usage:
#   print_array item1 item2 item3
# usage:
#   array=('item1' 'item2' 'item3')
#   print_array "${array[@]}"
# https://www.gnu.org/software/bash/manual/html_node/Arrays.html
print_array() {
  local list
  list=$(printf "'%s', " "${@}")
  echo "${list%, }"
}

# check that a command is available
# usage:
#   check_installed jq curl
check_installed() {
  # gather a list of missing commands
  missing_commands=()
  for cmd in "${@}"; do
    if ! command -v "$cmd" &>/dev/null; then
      missing_commands+=("$cmd")
    fi
  done

  # if any commands are missing, print an error and exit
  if [[ ${#missing_commands[@]} -gt 0 ]]; then
    local list
    list=$(print_array "${missing_commands[@]}")
    echo "Please install $list" >&2
    exit 1
  fi
}

# make a gitlab api call
# usage:
#   glab_api GET /projects/123/registry/repositories
#   glab_api POST /projects/123/registry/repositories "name=foo" "path=bar"
# needs: curl, $CI_API_V4_URL, $API_TOKEN
function glab_api {
  local method=${1:?}
  local endpoint=${2:?}
  local form_items=("${@:3}")

  local api_url=${CI_API_V4_URL:?}
  local api_token=${API_TOKEN:?}

  check_installed curl

  # construct the form data
  local form_data=""
  for item in "${form_items[@]}"; do
    form_data+="--form ${item} "
  done

  # curl options:
  # - no-progress-meter: don't show progress bar but show errors
  # - location: follow redirects
  # - fail: return 22 on server errors
  # - header: send private token for authentication
  # - request: specify the HTTP method
  # https://docs.gitlab.com/ee/api/container_registry.html
  echo "${api_url}${endpoint}"
  echo "${form_data}"
  curl \
    --no-progress-meter \
    --fail \
    --location \
    --header "PRIVATE-TOKEN: ${api_token}" \
    --request "${method}" \
    "${form_data}${api_url}${endpoint}"
}

# get the id of a container registry repo for a project by path
# usage:
#   get_container_registry_repo_id
# needs: jq, curl, $CI_API_V4_URL, $API_TOKEN, $CI_PROJECT_ID, $CI_PROJECT_PATH_SLUG
function get_container_registry_repo_id {
  local project_id=${CI_PROJECT_ID:?}
  local project_path_slug=${CI_PROJECT_PATH_SLUG:?}

  check_installed jq

  response=$(glab_api GET "/projects/${project_id}/registry/repositories")

  # jq query:
  # - raw-output: output raw strings instead of json
  # - exit-status: return 1 if no results
  # - for each elements of the array
  # - where the registry repo path is the one we care about
  # - select the id prop
  # https://jqlang.github.io/jq/manual/
  echo "${response}" | jq --raw-output --exit-status ".[] | select(.path == \"${project_path_slug}\") | .id"
}

# check if a tag exists in a project's container registry
# usage:
#   container_tag_exists 123 "my-tag"
# needs: curl, $CI_API_V4_URL, $API_TOKEN, $CI_PROJECT_ID
function container_tag_exists {
  local reg_repo_id=${1:?}
  local image_tag=${2:?}

  local project_id=${CI_PROJECT_ID:?}

  glab_api GET "/projects/${project_id}/registry/repositories/${reg_repo_id}/tags/${image_tag}" >/dev/null
}

# delete tag from a project's container registry
# usage:
#  delete_container_tag 123 "my-tag"
# needs: curl, $CI_API_V4_URL, $API_TOKEN, $CI_PROJECT_ID
function delete_container_tag {
  local reg_repo_id=${1:?}
  local image_tag=${2:?}

  local project_id=${CI_PROJECT_ID:?}

  glab_api DELETE "/projects/${project_id}/registry/repositories/${reg_repo_id}/tags/${image_tag}" >/dev/null
}

# build the container image with kaniko
# usage:
#   build_container "ubuntu:22.04"
# needs: /kaniko/executor, $HTTP_PROXY, $HTTPS_PROXY, $CI_REGISTRY_IMAGE
function build_container {
  local base_image=${1:?}

  local project_registry=${CI_REGISTRY_IMAGE:?}

  # strip tag from base image
  local base_image_no_tag=${base_image%%:*}

  # build with kaniko
  # https://github.com/GoogleContainerTools/kaniko
  /kaniko/executor \
    --context build \
    --build-arg "BASE_IMAGE=${base_image}" \
    --build-arg "http_proxy=${HTTP_PROXY}" \
    --build-arg "https_proxy=${HTTPS_PROXY}" \
    --dockerfile build/dockerfile \
    --destination "${project_registry}:${base_image_no_tag}-gridpack-env"
}

# build container image if not available in registry or force_rebuild = 'true'
# usage:
#   check_container_needs_built "ubuntu:22.04"
# needs: jq, curl, $CI_API_V4_URL, $API_TOKEN, $CI_PROJECT_ID, $CI_PROJECT_PATH_SLUG
function check_container_needs_built {
  local tag=${1:?}

  local repo_id
  repo_id=$(get_container_registry_repo_id)

  # if we see a file in the current directory called "container-definition-changed"
  # or the tag does not exist in the registry, then the container needs rebuilt
  if test -f ./container-definition-changed || ! container_tag_exists "$repo_id" "$tag"; then
    # in case the build fails, delete the image associated with this tag from the registry
    # this will signal to potential subsequent pipelines that the container image still needs
    # rebuilt even if relevant files have not changed
    # this provents unknowingly using an outdated image in subsequent pipelines
    delete_container_tag "$repo_id" "$tag"

    # signal to the next job via an empty file artifact that the container needs built
    touch ./container-needs-built
  fi
}
