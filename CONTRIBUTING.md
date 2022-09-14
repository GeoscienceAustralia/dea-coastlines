# Contributing to dea-coastlines

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to dea-coastlines and its packages. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

## Integration test data
### setting a NEW db dump

```
    # bring up indexing and db container
    docker-compose -f docker-compose.index.yaml -f docker-compose.cleandb.yaml up
```

### building on top of existing db dump

```
  docker-compose -f docker-compose.yaml -f docker-compose.index.yaml up
```

### indexing and create db dump

```
  # start by going to index container
  docker exec -ti dea-coastlines_index_1 bash
  datacube system init # OPTIONAL: no need to run this command if building off existing db
  ... # anything extra to index
  exit # after indexing is complete, exit

  # return to index container
  docker exec -it dea-coastlines_index_1 bash
  pg_dump -U localhost -p 5432 -h localhost odc > dump.sql
  # enter password on prompt: mysecretpassword or check .env file
  exit


  # copy the new dump to dea-coastlines/docker/database folder
  docker cp dea-coastlines_index_1:/dump.sql dea-coastlines/docker/database
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
