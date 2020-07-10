#!/usr/bin/env python

# Run with python setup.py develop --user

from setuptools import setup
import os
import platform

setup_requires = []
scripts = []

setup(
    name="popgen_utils",
    version="0.0.1",
    author='Lauren Sugden & Arthur Sugden',
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy>=1.11',
        'scipy>=0.18.0',
        'pandas>=0.18.0',
        'pyyaml',
        # 'scikit-learn>=0.11',
    ],
    scripts=scripts,
    # package_data={
    #     'replay': [
    #         'tests/*.py',
    #         'tests/data/example.tif',
    #         'tests/data/example.h5',
    #         'tests/data/example-volume.h5',
    #         'tests/data/example-tiffs/*.tif',
    #     ]
    # },
    #
    # metadata for upload to PyPI
    author_email="lauren.v.sugden@gmail.com",
    description="Sugden Lab population genetics analysis tools",
    license="GNU GPLv2",
    keywords="popgen evolutionary biology swifr",
    setup_requires=setup_requires,
    # setup_requires=['setuptools_cython'],
    platforms=["Linux", "Mac OS-X", "Windows"],
    #
    # could also include long_description, download_url, etc.
)
