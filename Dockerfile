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
    # For SSL
    ca-certificates \
    # Tidy up
    && apt-get autoclean && \
    apt-get autoremove && \
    rm -rf /var/lib/{apt,dpkg,cache,log}


# Environment can be whatever is supported by setup.py
# so, either deployment, test
ARG ENVIRONMENT=deployment
# ARG ENVIRONMENT=test

RUN echo "Environment is: $ENVIRONMENT"

COPY requirements.txt /tmp/
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r /tmp/requirements.txt \
    --no-binary rasterio \
    --no-binary shapely \
    --no-binary fiona \
    # Extras
    && pip install --no-cache-dir awscli requests

# Set up a nice workdir and add the live code
ENV APPDIR=/code
RUN mkdir -p $APPDIR
WORKDIR $APPDIR
ADD . $APPDIR

RUN if [ "$ENVIRONMENT" = "deployment" ] ; then\
        pip install .[$ENVIRONMENT] ; \
    else \
        pip install --editable .[$ENVIRONMENT] ; \
    fi


CMD ["python", "--version"]

RUN  deacoastlines-raster --help \
  && deacoastlines-vector --help