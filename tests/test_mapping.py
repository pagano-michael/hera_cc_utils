# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Define tests for mapping.py."""


import pytest
import numpy as np
from hera_cc_utils import mapping


def test_coords_in():
    """Define a test for coords_in."""
    assert mapping.coords_in("gsm") == "G"
    assert mapping.coords_in("roman") == "E"
    assert mapping.coords_in("ngrst") == "E"
    assert mapping.coords_in("euclid_south") == "C"
    assert mapping.coords_in("euclid_fornax") == "C"
    assert mapping.coords_in("so") == "C"

    return


def test_coords_in_warning():
    """Test that a warning is raised for an unknown dataset."""
    with pytest.warns(
        UserWarning, match="Assuming input map is in celestial coordinates"
    ):
        coords = mapping.coords_in("foo")

    assert coords == "C"

    return


def test_map_init_errors():
    """Test errors being raised in initializing Map object."""
    with pytest.raises(ValueError, match="Unrecognized input map: foo"):
        mapping.Map(data="foo")

    return


def test_coords_in_map():
    """Test the coords_in method on the Map object."""
    # try something for an existing data map
    mymap = mapping.Map(data="so")
    assert mymap.coords_in == "C"

    # try a stripe
    mymap = mapping.Map(data="hera")
    assert mymap.coords_in == "C"

    # try a user-defined thing
    data = np.arange(10)
    mymap = mapping.Map(data=data, coords_in="G")
    assert mymap.coords_in == "G"

    return


def test_is_binary():
    """Test the is_binary method on the Map object."""
    mymap = mapping.Map(data="so")
    assert mymap.is_binary

    data = np.arange(10)
    mymap = mapping.Map(data=data)
    assert not mymap.is_binary

    return
