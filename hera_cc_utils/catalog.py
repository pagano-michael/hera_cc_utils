# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Utilities for dealing with galaxy/QSO catalogs."""

import numpy as np
from astropy.coordinates import SkyCoord

from .util import deg_per_hr

_xshooter_ref = "https://ui.adsabs.harvard.edu/abs/2020ApJ...905...51S/abstract"

# VIKING
_viking_ref1 = "https://ui.adsabs.harvard.edu/abs/2013ApJ...779...24V/abstract"
_viking_ref2 = "https://ui.adsabs.harvard.edu/abs/2015MNRAS.453.2259V/abstract"
_viking = {
    "J2348–3054": {
        "ra": "23h48m33.34s",
        "dec": "-30d54m10.0s",
        "z": 6.886,
        "ref": _viking_ref1,
    },
    "J0109–3047": {
        "ra": "01h09m53.13s",
        "dec": "-30d47m26.3s",
        "z": 6.745,
        "ref": _viking_ref1,
    },
    "J0305–3150": {
        "ra": "03h05m16.92s",
        "dec": "-31d50m56.0s",
        "z": 6.604,
        "ref": _viking_ref1,
    },
    "J0328−3253": {
        "ra": "03h28m35.511s",
        "dec": "-32d53m22.92s",
        "z": 5.860,
        "ref": _viking_ref2,
    },
    "J0046–2837": {
        "ra": "00h46m23.645s",
        "dec": "-28d37m47.34s",
        "z": 5.9926,
        "ref": _xshooter_ref,
    },
    "J2211–3206": {
        "ra": "22h11m12.391s",
        "dec": "-32d06m12.95s",
        "z": 6.3394,
        "ref": _xshooter_ref,
    },
    "J2318–3029": {
        "ra": "23h18m33.103s",
        "dec": "-30d29m33.36s",
        "z": 6.1456,
        "ref": _xshooter_ref,
    },
    "J2348–3054_xshooter": {
        "ra": "23h48m33.336s",
        "dec": "-30d54m10.24s",
        "z": 6.9007,
        "ref": _xshooter_ref,
    },
}

# Pan-STARRS1
_ps1_ref1 = "https://ui.adsabs.harvard.edu/abs/2014AJ....148...14B/abstract"
_ps1_ref2 = "https://ui.adsabs.harvard.edu/abs/2017ApJ...849...91M/abstract"
_ps1 = {
    "PSO 231-20": {"ra": "231.6576", "dec": "-20.8335", "z": 6.5864, "ref": _ps1_ref2},
    "PSO J037.9706–28.8389": {
        "ra": "02h31m52.96s",
        "dec": "-28d50m20.1s",
        "z": 5.99,
        "ref": _ps1_ref1,
    },
    "PSO J065.4085–26.9543": {
        "ra": "04h21m38.049s",
        "dec": "-26d57m15.61s",
        "z": 6.1871,
        "ref": _xshooter_ref,
    },
}

# Banados+ 2016 https://ui.adsabs.harvard.edu/abs/2016ApJS..227...11B/abstract
# has table of all z > 5.6 quasars known at that point (March 2016).
# https://ned.ipac.caltech.edu/inrefcode?search_type=Search&refcode=2016ApJS..227...11B


# VLT ATLAS
# https://ui.adsabs.harvard.edu/abs/2015MNRAS.451L..16C/abstract
_atlas_ref1 = "https://ui.adsabs.harvard.edu/abs/2015MNRAS.451L..16C/abstract"
_atlas_ref2 = "https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.1649C/abstract"
_atlas = {
    "J025.6821-33.4627": {
        "ra": "025.6821",
        "dec": "-33.4627",
        "z": 6.31,
        "ref": _atlas_ref1,
    },
    "J332.8017−32.1036": {
        "ra": "332.8017",
        "dec": "-32.1036",
        "z": 6.32,
        "ref": _atlas_ref2,
    },
}

# VHS-DES
_ps1_vhs_des = "https://ui.adsabs.harvard.edu/abs/2019MNRAS.487.1874R/abstract"
_des = {
    "VDES J0020-3653": {
        "ra": "00h20m31.47s",
        "dec": "-36d53m41.8s",
        "z": 6.5864,
        "ref": _ps1_vhs_des,
    },
}

_yang = "https://ui.adsabs.harvard.edu/abs/2020ApJ...904...26Y/abstract"
_decarli = "https://ui.adsabs.harvard.edu/abs/2018ApJ...854...97D/abstract"
_other = {
    "J0142−3327": {"ra": "0142", "dec": "-3327", "z": 6.3379, "ref": _yang},
    "J0148−2826": {"ra": "0148", "dec": "-2826", "z": 6.54, "ref": _yang},
    "J2002−3013": {"ra": "2002", "dec": "-3013", "z": 6.67, "ref": _yang},
    "J2318–3113": {
        "ra": "23h18m18.351s",
        "dec": "−31d13m46.35s",
        "z": 6.444,
        "ref": _decarli,
    },
}


def _to_decimal(s):
    if "." in s:
        out = float(s)
    elif s[0] == "-":
        out = float(s[0:3] + "." + s[3:])
    else:
        out = float(s[0:2] + "." + s[2:])

    return out


_qso_catalogs = {"viking": _viking, "panstarrs": _ps1, "atlas": _atlas, "other": _other}


class Catalog(object):
    """
    Define a class for handling QSO catalogs.

    Parameters
    ----------
    data : str
        The type of data to handle. Right now "qso" is the only allowed value.
    kwargs : dict
        Keyword arguments to save directly on the object.
    """

    def __init__(self, data, **kwargs):
        self.data = data
        self.kwargs = kwargs

    def get_all_pos(self, zmin=None):
        """
        Return a list of (RA, DEC, redshift) for all objects.

        Parameters
        ----------
        zmin : float
            The minimum redshift to include for objects in the catalog.

        Returns
        -------
        names : list of str, shape (n_objects)
            The names of objects in the catalog.
        data : ndarray, shape (n_objects, 3)
            The RA [hour angle], dec [degree], and redshift of the objects.

        Raises
        ------
        ValueError
            This is raised if `self.data` is not "qso", as this is the only type
            of data we know how to handle right now.
        """
        if not self.data.lower().startswith("qso"):
            raise ValueError("Only know how to do QSOs right now.")

        data = []
        names = []
        for cat in _qso_catalogs.keys():

            for element in _qso_catalogs[cat]:
                obj = _qso_catalogs[cat][element]

                if zmin is not None:
                    if obj["z"] < zmin:
                        continue

                if "h" in obj["ra"]:
                    kw = {"frame": "icrs"}
                    ra = obj["ra"]
                    dec = obj["dec"]
                else:
                    kw = {"unit": "degree", "frame": "icrs"}
                    if len(obj["ra"]) == 4:
                        ra = _to_decimal(obj["ra"]) * deg_per_hr
                    else:
                        ra = _to_decimal(obj["ra"])

                    dec = _to_decimal(obj["dec"])

                coord = SkyCoord(ra, dec, **kw)

                names.append(element)
                data.append((coord.ra.hour, coord.dec.degree, obj["z"]))

        return names, np.array(data)
