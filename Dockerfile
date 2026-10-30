# Base image with:
# - Ubuntu 22.04
# - Python 3.10.12
# - GDAL 3.7.3, released 2023/10/30
FROM ghcr.io/osgeo/gdal:ubuntu-small-3.7.3

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Install system dependencies and remove the compromised
# emSign ECC Root CA - G3 certificate from the system trust store.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      openssl \
      build-essential \
      git \
      python3-dev \
      python3-pip \
      libpq-dev && \
    if grep -q '^mozilla/emSign_ECC_Root_CA_-_G3.crt$' /etc/ca-certificates.conf; then \
      sed -i \
        's|^mozilla/emSign_ECC_Root_CA_-_G3.crt$|!mozilla/emSign_ECC_Root_CA_-_G3.crt|' \
        /etc/ca-certificates.conf; \
    fi && \
    update-ca-certificates --fresh && \
    ! grep -Ril "emSign ECC Root CA - G3" \
      /etc/ssl/certs \
      /usr/local/share/ca-certificates \
      2>/dev/null && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# eMudhra/emSign CA compromise remediation.
ARG PIP_VERSION=26.2.2
ARG UV_VERSION=0.13.4

COPY requirements.in constraints.txt /app/

# Install pinned packaging tools.
RUN python3 -m pip install \
      --no-cache-dir \
      --upgrade \
      "pip==${PIP_VERSION}" && \
    python3 -m pip install \
      --no-cache-dir \
      "uv==${UV_VERSION}"

# Resolve dependencies subject to the security constraints.
RUN uv pip compile \
      --constraint /app/constraints.txt \
      /app/requirements.in \
      -o /app/requirements.txt

# Install resolved dependencies.
RUN uv pip install \
      --constraint /app/constraints.txt \
      -r /app/requirements.txt \
      --system

COPY . /app

# Install DEA Coastlines while continuing to enforce the
# eMudhra/emSign dependency constraints.
RUN uv pip install \
    . \
    --system \
    --constraint /app/constraints.txt
    
# Verify the resulting Python environment.
RUN uv pip check

# Verify DEA Coastlines entry points.
RUN deacoastlines-raster --help
RUN deacoastlines-vector --help
