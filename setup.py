"""
Setup script for RXNRECer
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Read requirements
requirements = []
with open("requirements.txt") as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith("#"):
            requirements.append(line)

setup(
    name="rxnrecer",
    version="1.0.0",
    author="Zhenkun Shi",
    author_email="kingstdio@gmail.com",
    description="Enzyme Reaction Prediction from Protein Sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/rxnrecer",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=22.0.0",
            "flake8>=5.0.0",
            "mypy>=1.0.0",
            "isort>=5.10.0",
            "pre-commit>=2.20.0",
        ],
        "docs": [
            "sphinx>=5.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "rxnrecer=rxnrecer.cli.predict:main",
        ],
    },
    include_package_data=True,
    package_data={
        "rxnrecer": ["config/*.json", "data/*"],
    },
    keywords="bioinformatics protein sequence enzyme reaction prediction deep learning",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/rxnrecer/issues",
        "Source": "https://github.com/yourusername/rxnrecer",
        "Documentation": "https://rxnrecer.readthedocs.io/",
    },
)
