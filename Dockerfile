# Base image with:
# - Ubuntu 22.04
# - Python 3.10.12
# - GDAL 3.7.3, released 2023/10/30
FROM ghcr.io/osgeo/gdal:ubuntu-small-3.7.3

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Apt installation
RUN apt-get update && \
    apt-get install -y \
      build-essential \
      fish \
      git \
      vim \
      htop \
      wget \
      unzip \
      python3-pip \
      libpq-dev \
    && apt-get autoclean && \
    apt-get autoremove && \
    rm -rf /var/lib/{apt,dpkg,cache,log}

# Install pip-tools
RUN pip install pip-tools

# Pip installation
RUN mkdir -p /conf
COPY requirements.txt /conf/
RUN pip install -r /conf/requirements.txt \
    && pip install --no-cache-dir awscli

# Copy source code and install it
RUN mkdir -p /code
WORKDIR /code
ADD . /code

RUN echo "Installing dea-coastlines through the Dockerfile."
RUN pip install --extra-index-url="https://packages.dea.ga.gov.au" .

RUN pip freeze && pip check

# Make sure it's working
RUN deacoastlines-raster --help \
  && deacoastlines-vector --help
