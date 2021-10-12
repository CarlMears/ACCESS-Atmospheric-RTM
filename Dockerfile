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
# docker run access-atmospheric-rtm:latest -m access_atmosphere.download ...
#
# Or:
#
# docker run access-atmospheric-rtm:latest -m access_atmosphere.process...
#
# Probably a bind-mount volume is needed to store the outputs. For example:
#
# docker run -v $PWD:/data access-atmospheric-rtm:latest -m access_atmosphere.process /data/era5_surface.nc /data/era5_levels.nc /data/rtm_out.nc

FROM docker.io/library/python:3.9 AS build

RUN apt-get update && \
    apt-get install -y gfortran
RUN python3 -m venv --upgrade-deps /root/venv
RUN /root/venv/bin/pip3 install build

WORKDIR /root
COPY . .

RUN /root/venv/bin/python3 -m build

FROM docker.io/library/python:3.9-slim

RUN apt-get update && \
    apt-get install -y libgfortran5 libgomp1
RUN python3 -m venv --upgrade-deps /root/venv

WORKDIR /root
COPY --from=build \
    /root/dist/*.whl \
    /root/

RUN /root/venv/bin/pip3 install ./*.whl
ENTRYPOINT [ "/root/venv/bin/python3" ]

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
