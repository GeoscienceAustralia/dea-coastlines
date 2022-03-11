from logging import Logger, basicConfig, INFO
import yaml
import fsspec


def configure_logging(name: str = "Coastlines") -> Logger:
    """
    Configure logging for the application.
    """
    basicConfig(
        format="%(asctime)s %(levelname)s %(message)s",
        level=INFO,
        filename="{}.log".format(name),
        filemode="w",
    )
    # Create a logger and return it
    return Logger(name)


def load_config(config_path: str) -> dict:
    """
    Loads a YAML config file and returns data as a nested dictionary.

    config_path can be a path or URL to a web accessible YAML file
    """
    with fsspec.open(config_path, mode="r") as f:
        config = yaml.safe_load(f)
    return config
