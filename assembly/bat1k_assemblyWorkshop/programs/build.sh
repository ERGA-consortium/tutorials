#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'

VERSION=v0.6.3
DOCKER_TAG="assembly-workshop:$VERSION"
SINGULARITY_IMAGE="${DOCKER_TAG/:/_}.sif"

# build Docker image
docker build -t "$DOCKER_TAG" -f Dockerfile .

# test Docker image
#$docker run --rm "$DOCKER_TAG"

# build Singularity image
SINGULARITY_DISABLE_CACHE=true singularity build "$SINGULARITY_IMAGE" "docker-daemon://$DOCKER_TAG"

# test Singularity image
$singularity exec --no-home --cleanenv assembly-workshop_${VERSION}.sif /bin/bash 
