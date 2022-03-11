import logging
import yaml
import fsspec


def configure_logging(name: str = "Coastlines") -> logging.Logger:
    """
    Configure logging for the application.
    """
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO
    )
    # Create a logger and return it
    return logging.getLogger(name)


def load_config(config_path: str) -> dict:
    """
    Loads a YAML config file and returns data as a nested dictionary.

    config_path can be a path or URL to a web accessible YAML file
    """
    with fsspec.open(config_path, mode="r") as f:
        config = yaml.safe_load(f)
    return config
