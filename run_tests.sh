#!/bin/bash

# Can never remember all the flags
pytest --cov-config=.coveragerc --cov=hera_cc_utils --cov-report=html -v tests/*.py
