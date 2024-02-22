The shell scripts in this directory are used to automate the installation of GridPACK on Debian/RHEL based systems.

# GitLab Mirror

A [PNNL GitLab instance](https://devops.pnnl.gov/gridpack-code/GridPACK) hosting a [pull-mirror](https://docs.gitlab.com/ee/user/project/repository/mirror/pull.html) of this [GitHub repo](https://github.com/GridOPTICS/GridPACK) uses these scripts to execute tests whenever commits are pushed to GitHub. GitLab polls the GitHub repo to check for new commits every 30 minutes or every 5 minutes when triggered by the GitLab UI/API. An [integration](https://docs.gitlab.com/ee/ci/ci_cd_for_external_repos/github_integration.html) was configured to report the result back to GitHub for any commit associated with a build ["pipeline"](https://docs.gitlab.com/ee/ci/pipelines/).

# GitLab Jobs

The test automation is defined by the following set of ["jobs"](https://docs.gitlab.com/ee/ci/jobs/) in `/.gitlab-ci.yml`, which may be added to a pipeline by GitLab when a build is triggered:

1. `check-container-definition-changed`: Creates a file called `container-definition-changed` if the scripts or `dockerfile` used to build the container image have been modified. The file is passed to subsequent jobs as an ["artifact"](https://docs.gitlab.com/ee/ci/jobs/job_artifacts.html).
2. `.check-container-needs-built`: A template job for multiple distributions. The leading `.` on this job name tells GitLab it is not meant to be added to the pipeline. This job runs `check_container_needs_built`, a function in `container_lib.sh` which creates another signal file called `container-needs-built` if the definition has changed or the container image tag does not exist in the GitLab project [container registry](https://docs.gitlab.com/ee/user/packages/container_registry/). Additionally, if the container image does need to be rebuilt, this function deletes the image tag from the registry to avoid unintentional use by subsequent pipelines in the event that the container image build fails to complete successfully and overwrite the old image.
3. `.build-container`: A template job for multiple distributions. Uses [kaniko](https://github.com/GoogleContainerTools/kaniko) to build a container image based on a given distribution if the file `container-needs-built` is passed as an artifact from a prior job or an environment variable called `FORCE_CONTAINER_BUILD` is set to `true`. The environment variable is intended to be used when [manually triggering a pipeline from the GitLab UI](https://docs.gitlab.com/ee/ci/pipelines/#run-a-pipeline-manually).
4. `.build-gridpack`: A template job for multiple distributions. Builds GridPACK in a container using an image from the GitLab project container registry based on a given distribution whenever files in `/src`, `/python`, or `/build` are modified using `install_gridpack.sh`. All files created by the build are saved as artifacts.
5. `.test-gridpack`: A template job for multiple distributions. Tests GridPACK in a new container using the  artifacts from the associated `.build-container` job and with the same container image. These tests could have been performed at the end of the `.build-container` job but were separated to provide a distinct status for each concern and to maintain organized logs/artifacts.

# Container image

`dockerfile` defines a container image used to create an environment with all necessary GridPACK dependencies, including those installed from packages and those installed from source.
