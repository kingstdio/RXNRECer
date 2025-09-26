"""
Version information for RXNRECer.
This file serves as the single source of truth for version numbers.
"""

__version__ = "1.3.5"
__version_info__ = (1, 3, 5)

# Version components
MAJOR_VERSION = 1
MINOR_VERSION = 3
PATCH_VERSION = 0

# Build information
BUILD_DATE = "2025-09-26"
BUILD_COMMIT = "53b5ffe"

# Release information
RELEASE_NAME = "Enhanced Features and Optimizations"
RELEASE_DATE = "2024-09-26"

def get_version():
    """Get the version string."""
    return __version__

def get_version_info():
    """Get the version tuple."""
    return __version_info__

def get_full_version():
    """Get the full version string with build info."""
    return f"{__version__} ({BUILD_DATE})"

def get_release_info():
    """Get release information."""
    return {
        "version": __version__,
        "version_info": __version_info__,
        "build_date": BUILD_DATE,
        "build_commit": BUILD_COMMIT,
        "release_name": RELEASE_NAME,
        "release_date": RELEASE_DATE
    }
