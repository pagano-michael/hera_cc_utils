# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Tests for line.py."""

import pytest
import numpy as np
from astropy import constants
from astropy import units

from hera_cc_utils import line
from hera_cc_utils.util import line_registry


@pytest.mark.parametrize("line_keyvals", line_registry.items())
def test_get_line_name(line_keyvals):
    """Test get_line_name method."""
    # get line from line_registry
    line_name, line_dict = line_keyvals
    if "freq_rest" in line_dict:
        li = line.Line(line_name, freq_rest=line_dict["freq_rest"])
    else:
        li = line.Line(line_name, lambda_rest=line_dict["lambda_rest"])

    assert line_name == li.get_line_name()

    return


@pytest.mark.parametrize("line_keyvals", line_registry.items())
def test_get_rest_lambda(line_keyvals):
    """Test get_rest_lambda method."""
    line_name, line_dict = line_keyvals
    if "freq_rest" in line_dict:
        li = line.Line(line_name, freq_rest=line_dict["freq_rest"])
        lambda_rest = (
            (constants.c / (line_dict["freq_rest"] * units.MHz)).to("um").value
        )  # convert to µm and extract numerical value
    else:
        li = line.Line(line_name, lambda_rest=line_dict["lambda_rest"])
        lambda_rest = line_dict["lambda_rest"]

    assert np.isclose(li.get_rest_lambda(), lambda_rest)

    return


@pytest.mark.parametrize("line_keyvals", line_registry.items())
def test_get_rest_freq(line_keyvals):
    """Test get_rest_freq method."""
    line_name, line_dict = line_keyvals
    if "freq_rest" in line_dict:
        li = line.Line(line_name, freq_rest=line_dict["freq_rest"])
        freq_rest = line_dict["freq_rest"]
    else:
        li = line.Line(line_name, lambda_rest=line_dict["lambda_rest"])
        freq_rest = (
            (constants.c / (line_dict["lambda_rest"] * units.um)).to("MHz").value
        )  # convert to MHz and extract numerical value

    assert np.isclose(li.get_rest_freq(), freq_rest)

    return


@pytest.mark.parametrize("line_keyvals", line_registry.items())
def test_z_to_lambda(line_keyvals):
    """Test z_to_lambda method."""
    line_name, line_dict = line_keyvals
    if "freq_rest" in line_dict:
        li = line.Line(line_name, freq_rest=line_dict["freq_rest"])
    else:
        li = line.Line(line_name, lambda_rest=line_dict["lambda_rest"])

    lambda_rest = li.get_rest_lambda()

    zvals = np.arange(20)
    for z in zvals:
        reference_lambda = lambda_rest * (1 + z)
        assert np.isclose(reference_lambda, li.z_to_lambda(z))

    return


@pytest.mark.parametrize("line_keyvals", line_registry.items())
def test_z_to_freq(line_keyvals):
    """Test z_to_freq method."""
    line_name, line_dict = line_keyvals
    if "freq_rest" in line_dict:
        li = line.Line(line_name, freq_rest=line_dict["freq_rest"])
    else:
        li = line.Line(line_name, lambda_rest=line_dict["lambda_rest"])

    freq_rest = li.get_rest_freq()

    zvals = np.arange(20)
    for z in zvals:
        reference_freq = freq_rest / (1 + z)
        assert np.isclose(reference_freq, li.z_to_freq(z))

    return


@pytest.mark.parametrize("line_keyvals", line_registry.items())
def test_lambda_to_z(line_keyvals):
    """Test lambda_to_z method."""
    line_name, line_dict = line_keyvals
    if "freq_rest" in line_dict:
        li = line.Line(line_name, freq_rest=line_dict["freq_rest"])
    else:
        li = line.Line(line_name, lambda_rest=line_dict["lambda_rest"])

    lambda_rest = li.get_rest_lambda()

    lambda_obs = lambda_rest * 10
    assert np.isclose(li.lambda_to_z(lambda_obs), 9.0)

    return


@pytest.mark.parametrize("line_keyvals", line_registry.items())
def test_freq_to_z(line_keyvals):
    """Test freq_to_z method."""
    line_name, line_dict = line_keyvals
    if "freq_rest" in line_dict:
        li = line.Line(line_name, freq_rest=line_dict["freq_rest"])
    else:
        li = line.Line(line_name, lambda_rest=line_dict["lambda_rest"])

    freq_rest = li.get_rest_freq()

    freq_obs = freq_rest / 10
    assert np.isclose(li.freq_to_z(freq_obs), 9.0)

    return
