# Integration test setup

## Integration test data
### Setting up a new database dump

```
    # Bring up indexing and database container
    docker-compose -f docker-compose.index.yaml -f docker-compose.cleandb.yaml up

    # Start by going to index container
    docker exec -ti dea-coastlines_index_1 bash
    datacube system init
    exit
```

### Building on top of existing database dump

```
  docker-compose -f docker-compose.yaml -f docker-compose.index.yaml up
```

### Indexing and creating database dump

The products indexed into the existing database dump are:

- https://explorer.dev.dea.ga.gov.au/product/ga_ls8c_ard_3/regions/089083
- https://explorer.dev.dea.ga.gov.au/product/ga_ls5t_ard_3/regions/089083
- https://explorer.dev.dea.ga.gov.au/product/ga_ls7e_ard_3/regions/089083
- https://explorer.dev.dea.ga.gov.au/product/ga_ls9c_ard_3/regions/089083

```
  # Start by going to index container
  docker exec -ti dea-coastlines_index_1 bash

  # Indexing example:
  # - Add a new product
  datacube product add https://raw.githubusercontent.com/GeoscienceAustralia/dea-config/master/products/baseline_satellite_data/c3/ga_ls9c_ard_3.odc-product.yaml

  # - Index datasets from s3
  s3-to-dc --skip-lineage s3://dea-public-data/baseline/ga_ls9c_ard_3/089/083/**/*.odc-metadata.yaml ga_ls9c_ard_3

  pg_dump -U odc -p 5432 -h postgres odc > dump.sql
  # Enter password on prompt: odcpass or to check echo $DB_PASSWORD

  # Copy the new database dump to the dea-coastlines/docker/database folder
  docker cp dea-coastlines_index_1:/dump.sql dea-coastlines/docker/db/dump.sql
```

### Running tests locally

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
