#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*-

from setuptools import setup, find_namespace_packages

link = "https://github.com/HERA-Team/hera_cc_utils"

setup_args = {
    "name": "hera_cc_utils",
    "version": "0.1",
    "description": "utilities for work on cross-correlations",
    "author": "Jordan Mirocha",
    "author_email": "mirochaj@gmail.com",
    "url": link,
    "packages": find_namespace_packages(),
    "include_package_data": True,
    "use_scm_version": True,
    "install_requires": [
        "healpy",
        "numpy",
        "matplotlib",
        "setuptools_scm",
    ],
}


setup(**setup_args)
