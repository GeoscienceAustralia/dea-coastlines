# Contributing to dea-coastlines

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to dea-coastlines and its packages. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

## Integration test data
### setting a NEW db dump

```
    # bring up indexing and db container
    docker-compose -f docker-compose.index.yaml -f docker-compose.cleandb.yaml up

    # start by going to index container
    docker exec -ti dea-coastlines_index_1 bash
    datacube system init
    exit
```

### building on top of existing db dump

```
  docker-compose -f docker-compose.yaml -f docker-compose.index.yaml up
```

### indexing and create db dump

The products indexed in existing db dump are:

- https://explorer.dev.dea.ga.gov.au/product/ga_ls8c_ard_3/regions/089083,
- https://explorer.dev.dea.ga.gov.au/product/ga_ls5t_ard_3/regions/089083,
- https://explorer.dev.dea.ga.gov.au/product/ga_ls7e_ard_3/regions/089083
- https://explorer.dev.dea.ga.gov.au/product/ga_ls9c_ard_3/regions/089083

```
  # start by going to index container
  docker exec -ti dea-coastlines_index_1 bash

  # indexing example
  # add a new product
  datacube product add https://raw.githubusercontent.com/GeoscienceAustralia/dea-config/master/products/baseline_satellite_data/c3/ga_ls9c_ard_3.odc-product.yaml

  # index datasets from s3
  s3-to-dc --skip-lineage s3://dea-public-data/baseline/ga_ls9c_ard_3/089/083/**/*.odc-metadata.yaml ga_ls9c_ard_3

  pg_dump -U odc -p 5432 -h postgres odc > dump.sql
  # enter password on prompt: odcpass or to check echo $DB_PASSWORD

  # copy the new dump to dea-coastlines/docker/database folder
  docker cp dea-coastlines_index_1:/dump.sql dea-coastlines/docker/db/dump.sql
```

### Running test locally

```
    mkdir artifacts
    chmod a+rw artifacts
    wget https://www.dropbox.com/s/ivx93rcdl9yfdaf/tide_models_clipped.zip?dl=1 -O tide_models_clipped.zip
    unzip tide_models_clipped.zip

    docker-compose up -d
    docker-compose exec -T coastline /bin/sh -c "sh ./docker/coastline/wait-for-db; pytest --cov=dea_coastlines --cov-report=xml tests/"
    docker-compose exec -T coastline /bin/sh -c "cp /tmp/coverage.xml /mnt/artifacts"
    docker-compose down
```
