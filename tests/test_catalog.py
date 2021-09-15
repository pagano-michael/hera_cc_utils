# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Define tests for catalog.py."""

import pytest
import numpy as np
from hera_cc_utils import catalog


def test_to_decimal():
    """Define a test for _to_decimal."""
    mystr = "37.5"
    assert np.isclose(catalog._to_decimal(mystr), 37.5)

    # test passing in 4-digit coordinates
    mystr = "0142"
    assert np.isclose(catalog._to_decimal(mystr), 1.42)

    # now a negative coordinate
    mystr = "-3327"
    assert np.isclose(catalog._to_decimal(mystr), -33.27)

    return


def test_get_all_pos():
    """Define a test for get_all_pos."""
    indata = catalog._qso_catalogs
    cat = catalog.Catalog("qso")
    names, outdata = cat.get_all_pos()

    for sourcename, data in indata.items():
        for obj in data.keys():
            assert obj in names

    return


def test_get_all_pos_zmin():
    """Define a test where we impose a minimum redshift cut."""
    cat = catalog.Catalog("qso")
    zcuts = [6, 6.5, 7]
    # Note! These numbers will change if we add more objects to the catalog
    nobjects = [14, 7, 0]
    for zmin, nobj in zip(zcuts, nobjects):
        names, outdata = cat.get_all_pos(zmin=zmin)

        assert len(names) == nobj
        # make sure all of our redshifts are >= minimum value
        if nobj > 0:
            assert not np.any(outdata[:, 2] < zmin)

    return


def test_get_all_pos_error():
    """Check that we get an error when we don't specify qso."""
    cat = catalog.Catalog("foo")
    with pytest.raises(ValueError, match="Only know how to do QSOs right now"):
        cat.get_all_pos()

    return
