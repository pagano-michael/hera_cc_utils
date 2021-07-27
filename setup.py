#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*-

import os
import subprocess
from setuptools import setup, Extension

link = "https://github.com/HERA-Team/hera_cc_utils"

setup_args = {
    "name": "hera_cc_utils",
    "version": "0.1",
    "description": "utilities for work on cross-correlations",
    "author": "Jordan Mirocha",
    "author_email": "mirochaj@gmail.com",
    "url": link,
    "packages": ["hera_cc_utils"],
}


setup(**setup_args)
