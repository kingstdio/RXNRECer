"""
Command line interface module for RXNRECer.
"""

def main(*args, **kwargs):
    from .predict import main as _main

    return _main(*args, **kwargs)


def download_data(*args, **kwargs):
    from .download import download_data as _download_data

    return _download_data(*args, **kwargs)


def cache(*args, **kwargs):
    from .cache import cache as _cache

    return _cache(*args, **kwargs)

__all__ = ["main", "download_data", "cache"]
