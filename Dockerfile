# ACCESS atmospheric RTM Dockerfile
#
# Build with:
#
# docker build -t access-atmospheric-rtm:latest \
#   --build-arg=build_date=$(date --utc --date="@${SOURCE_DATE_EPOCH:-$(date +%s)}" --rfc-3339='seconds') \
#   --build-arg=source_date=$(git log -1 --pretty=%ci) \
#   --build-arg=version=$(git rev-parse --short HEAD) \
#   --build-arg=revision=$(git symbolic-ref --short HEAD) \
#   .
#
# To run it:
#
# TODO: finish the running instructions. Note that the only bind-mount needed is
# /mnt/ops1p/n for NCEP data. Plus either a bind-mount for the output data or a
# volume.

FROM docker.io/library/python:3.9-slim AS build

RUN apt-get update && \
    apt-get install -y gfortran
RUN pip install --upgrade pip && pip install build

WORKDIR /root
COPY . .

RUN python -m build

FROM docker.io/library/python:3.9-slim

RUN pip install --upgrade pip

WORKDIR /root
COPY --from=build \
    /root/dist/*.whl \
    /root/

RUN pip install ./*.whl

ARG version
ARG revision
ARG build_date
ARG source_date
LABEL Maintainer="Richard Lindsley <lindsley@remss.com>" \
    "org.opencontainers.image.title"="ACCESS Atmospheric RTM" \
    "org.opencontainers.image.description"="Microwave atmospheric RTM for ACCESS work" \
    "org.opencontainers.image.created"="$build_date" \
    "org.opencontainers.image.authors"="Richard Lindsley <lindsley@remss.com>" \
    "org.opencontainers.image.vendor"="Remote Sensing Systems" \
    "org.opencontainers.image.url"="http://gitlab.remss.com/access/atmospheric-rtm" \
    "org.opencontainers.image.source"="http://gitlab.remss.com/access/atmospheric-rtm" \
    "org.opencontainers.image.version"="$version" \
    "org.opencontainers.image.revision"="$revision"
