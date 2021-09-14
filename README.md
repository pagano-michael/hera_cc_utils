# hera_cc_utils

A collection of modules to help determine the overlap between different surveys,
intensity mapping and otherwise, in real and Fourier space.

![Tests](https://github.com/HERA-Team/hera_cc_utils/actions/workflows/test_suite.yaml/badge.svg) [![codecov](https://codecov.io/gh/HERA-Team/hera_cc_utils/branch/main/graph/badge.svg?token=18ZMZEUWPW)](https://codecov.io/gh/HERA-Team/hera_cc_utils)

# Installation

After cloning the repository, installation can be handled by `pip` in the
typical fashion:

```bash
pip install .
```

## Dependencies

The following packages are required:

* astropy
* numpy
* matplotlib
* setuptools_scm

The following packages are optional:

* cartopy
* [pygsm](https://github.com/telegraphic/PyGSM)

# Testing

The repo uses `pytest` for the test suite. Simply run:

```bash
pytest
```
