FROM osgeo/gdal:ubuntu-small-3.4.1 as base

ENV CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt

RUN apt-get update \
    && apt-get install -y \
    # Build tools
    build-essential \
    git \
    python3-pip \
    # For Psycopg2
    libpq-dev python3-dev \
    # For mpi4py
#     libopenmpi-dev \
    # For cartopy
#     libproj-dev \
#     libgeos-dev \
    # For SSL
    ca-certificates \
    # Tidy up
    && apt-get autoclean && \
    apt-get autoremove && \
    rm -rf /var/lib/{apt,dpkg,cache,log}

COPY requirements.txt /tmp/
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r /tmp/requirements.txt \
    --no-binary rasterio \
    --no-binary shapely \
    --no-binary fiona \
    # Extras
    && pip install --no-cache-dir awscli requests


RUN mkdir -p /code
WORKDIR /code

COPY . /code/

RUN pip install /code

CMD ["python", "--version"]

RUN  deafricacoastlines-raster --help \
  && deafricacoastlines-vector --help