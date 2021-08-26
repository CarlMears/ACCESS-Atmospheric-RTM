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

FROM quay.io/rockylinux/rockylinux:8 AS build

RUN dnf install -y dnf-plugins-core epel-release gcc-gfortran gcc-c++ && \
    dnf config-manager --set-enabled powertools && \
    dnf install -y meson ninja-build netcdf-fortran-devel

WORKDIR /root
COPY rtm rtm

RUN meson setup --buildtype release build/ rtm/
RUN meson compile -C build/

FROM quay.io/rockylinux/rockylinux:8

RUN dnf install -y dnf-plugins-core epel-release && \
    dnf config-manager --set-enabled powertools && \
    dnf install -y --nodocs libgfortran libgomp netcdf-fortran python39 python39-pip python39-wheel && \
    dnf clean all && \
    rm -rf /var/cache/dnf/*

RUN pip3 install cdsapi

WORKDIR /root

COPY --from=build \
    /root/build/access_rtm \
    /root/build/access_rtm_ncep \
    /usr/local/bin/

COPY download_era5.py .

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
    "org.opencontainers.image.url"="http://gitlab.remss.com/lindsley/access-atmospheric-rtm" \
    "org.opencontainers.image.source"="http://gitlab.remss.com/lindsley/access-atmospheric-rtm" \
    "org.opencontainers.image.version"="$version" \
    "org.opencontainers.image.revision"="$revision"
