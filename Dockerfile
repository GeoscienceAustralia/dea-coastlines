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
      git \
      python3-pip \
      libpq-dev \
    && apt-get autoclean && \
    apt-get autoremove && \
    rm -rf /var/lib/{apt,dpkg,cache,log}

# Set up working directory
WORKDIR /app

# Copy requirements file first
COPY requirements.in /app/requirements.in

# Install uv and requirements
RUN pip install uv && \
    uv pip compile /app/requirements.in -o /app/requirements.txt && \
    uv pip install -r /app/requirements.txt --system

# Now copy the rest of the files
COPY . /app

# Install DEA Intertidal and verify installation
RUN uv pip install . --system && \
    uv pip check && \
    deacoastlines-raster --help && \
    deacoastlines-vector --help
