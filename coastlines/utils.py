import logging
import yaml
import fsspec
from pathlib import Path

STYLES_FILE = Path(__file__).parent / "styles.csv"


def configure_logging(name: str = "Coastlines") -> logging.Logger:
    """
    Configure logging for the application.
    """
    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    return logger


def load_config(config_path: str) -> dict:
    """
    Loads a YAML config file and returns data as a nested dictionary.

    config_path can be a path or URL to a web accessible YAML file
    """
    with fsspec.open(config_path, mode="r") as f:
        config = yaml.safe_load(f)
    return config
