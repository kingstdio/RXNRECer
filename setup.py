#!/usr/bin/env python3
"""
Setup script for RXNRECer package.
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="rxnrecer",
    version="1.0.0",
    author="Zhenkun Shi",
    author_email="zhenkun.shi@tib.cas.cn",
    description="Deep learning framework for predicting enzyme-catalyzed reactions from protein sequences",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/kingstdio/RXNRECer",
    project_urls={
        "Bug Reports": "https://github.com/kingstdio/RXNRECer/issues",
        "Source": "https://github.com/kingstdio/RXNRECer",
        "Documentation": "https://github.com/kingstdio/RXNRECer#readme",
    },
    packages=find_packages(include=['rxnrecer', 'rxnrecer.*']),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=22.0.0",
            "flake8>=5.0.0",
            "isort>=5.10.0",
        ],
        "docs": [
            "sphinx>=5.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "rxnrecer=rxnrecer.cli.predict:main",
            "rxnrecer-download-data=rxnrecer.cli.download:download_data",
            "rxnrecer-cache=rxnrecer.cli.cache:cache",
        ],
    },
    include_package_data=True,
    package_data={
        "rxnrecer": ["config/*.yaml", "config/*.yml"],
    },
    zip_safe=False,
    keywords="bioinformatics, machine-learning, deep-learning, protein, enzyme, reaction-prediction",
)
