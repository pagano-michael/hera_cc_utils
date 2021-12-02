# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Utilities for dealing with maps."""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import healpy

from .util import field_to_healpix
from .data import PATH


map_options = ["gsm", "roman", "ngrst", "euclid_south", "euclid_fornax", "so"]

mollview_projections = ["galactic", "equatorial", "celestial", "ecliptic"]

_coords_in = {
    "gsm": "G",
    "roman": "E",
    "ngrst": "E",
    "euclid_south": "C",
    "euclid_fornax": "C",
    "so": "C",
}

_filenames = {
    "roman": "Roman_GRISM_mask.fits",
    "euclid_south": "Euclid_deep_south.fits",
    "euclid_fornax": "Euclid_deep_Fornax.fits",
    "so": "so_lat.fits",
}

_stripes = {
    "hera": (-35, -25, 0, 24),
    # From Aihara et al. 2018, Table 5.
    "hsc_north": (42, 44.5, "13h20m00s", "16h40m00s"),
    "hsc_spring": (-2, 5, "09h00m00s", "15h30m00s"),
    "hsc_fall": (-1, 7, "22h00m00s", "02h40m00s"),
}


def coords_in(data):
    """
    Define the input coordinates to use for the given data.

    Parameters
    ----------
    data : str
        The dataset to use. Should be a key corresponding to the _coords_in
        dict.

    Returns
    -------
    str
        The coordinate system for the data. Should be one of: "G" (galactic),
        "E" (ecliptic), or "C" (celestial or equatorial).
    """
    if isinstance(data, str):
        if data in _coords_in:
            return _coords_in[data]

    warnings.warn(
        "Assuming input map is in celestial coordinates. Set `coords_in` by "
        "hand to 'G', 'C', or 'E' to change."
    )
    return "C"


class Map(object):
    """
    Utilities for dealing with maps.

    Right now, these utilities assume all-sky maps like the global sky model
    (GSM).

    .. note :: Can also be used to represent survey footprints, i.e., healpix
        maps with just ones and zeros.

    Parameters
    ----------
    data : str
        Data source for map. Current options include: 'gsm'.
    freq_unit : str
        Assumed frequency units when generating maps. Options: Hz, MHz, GHz.
    kwargs : dict
        Any other optional keyword arguments accepted by the map being
        generated, e.g., that can be supplied to pygsm.GlobalSkyModel.
    """

    def __init__(self, data="gsm", freq_unit="Hz", coords_in=None, nside=512, **kwargs):

        if isinstance(data, str):
            if (data not in map_options) and (data not in _stripes.keys()):
                raise ValueError(f"Unrecognized input map: {data}")

        if isinstance(data, str) and data in _stripes.keys():
            self.data = field_to_healpix(*_stripes[data], nside=nside)
            self._coords_in_ = "C"
        else:
            self.data = data
            self._coords_in_ = coords_in

        self.freq_unit = freq_unit
        self.kwargs = kwargs

    @property
    def coords_in(self):
        """
        Get the coordinates for the map.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The coordinate system for the map. One of: "G" (galactic), "E"
            (ecliptic), "C" (celestial or equatorial).
        """
        if not hasattr(self, "_coords_in"):
            if self._coords_in_ is None:
                self._coords_in = coords_in(self.data)
            else:
                self._coords_in = self._coords_in_

        return self._coords_in

    @property
    def is_binary(self):
        """
        Define whether the map is binary or not.

        Parameters
        ----------
        None

        Returns
        -------
        bool
           Whether the map is binary or not.
        """
        if not hasattr(self, "_is_binary"):
            if type(self.data) == np.ndarray:
                self._is_binary = np.unique(self.data).size == 2
            else:
                self._is_binary = self.data != "gsm"

        return self._is_binary

    @property
    def _gsm(self):
        """Return an instance of pygsm.GlobalSkyModel."""
        try:
            from pygsm import GlobalSkyModel
        except ImportError:
            raise ImportError(
                "Must have pygsm installed for working with the global sky model"
            )

        if not hasattr(self, "_gsm_"):
            self._gsm_ = GlobalSkyModel(freq_unit=self.freq_unit, **self.kwargs)
        return self._gsm_

    def get_map(self, freq=150e6, projection="ecliptic"):
        """
        Get a map suitable for plotting.

        .. note :: Currently this isn't general at all, just designed to work
            with GSM. Could be generalized to, e.g., plot Planck maps in the
            background.

        Parameters
        ----------
        freq : float, optional
            The frequency of the map, in Hz.
        projection : str, optional
            The projection to use. Must be one of: "galactic", "ecliptic",
            "celestial", "equatorial", "rectilinear", "Robinson"

        Returns
        -------
        img : ndarray
            The HEALPix map corresponding to the joint map.

        Raises
        ------
        ValueError
            Raised if `projection` is not valid.
        """
        coord_in = self.coords_in

        rot_kw = {"deg": False, "inv": True}
        if projection == "galactic":
            coord_out = "G"
        elif projection == "ecliptic":
            coord_out = "E"
        elif projection in ["celestial", "equatorial"]:
            coord_out = "C"
        elif projection in ["rectilinear", "Robinson"]:
            coord_out = "C"
        else:
            raise ValueError("Unrecognized option for `projection`.")

        if type(self.data) == np.ndarray:
            raw_map = self.data
            nside = healpy.npix2nside(raw_map.size)
        elif self.data in ["gsm"]:
            raw_map = self._gsm.generate(freq)
            nside = 512
        elif self.data in ["roman", "ngrst"]:
            fn = _filenames["roman"]
            raw_map = healpy.fitsfunc.read_map("{}/roman/{}".format(PATH, fn))
            nside = healpy.npix2nside(raw_map.size)
        elif self.data in ["euclid_south", "euclid_fornax"]:
            fn = _filenames[self.data]
            raw_map = healpy.fitsfunc.read_map("{}/euclid/{}".format(PATH, fn))
            nside = healpy.npix2nside(raw_map.size)
        elif self.data == "so":
            fn = _filenames[self.data]
            raw_map = healpy.fitsfunc.read_map("{}/so/{}".format(PATH, fn))
            nside = healpy.npix2nside(raw_map.size)
        elif type(self.data) == np.ndarray:
            raw_map = self.data
            nside = healpy.npix2nside(raw_map.size)
        else:
            # Map is not of type np.ndarray 
            raise NotImplementedError("Map is not of type np.ndarray or is not supported")

        # Initialize Rotator object to transform coordinates.
        rot = healpy.Rotator(coord=[coord_in, coord_out], **rot_kw)

        # Get (theta, phi) angles for all pixels
        theta, phi = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)))

        # Rotate into desired coordinate system.
        if rot is not None:
            theta_cel, phi_cel = rot(theta, phi)
            map_theta_phi = healpy.get_interp_val(raw_map, theta_cel, phi_cel)
        else:
            theta_cel, phi_cel = theta, phi
            map_theta_phi = raw_map

        # setup RA Dec coords

        # interpolate healpix
        if projection == "rectilinear":
            # Interpolate onto a rectangular mesh with size (nside, nside)
            # ra = np.linspace(-90, 270, nside)
            ra = np.linspace(0, 360, nside)
            dec = np.linspace(0, 180, nside)

            ra, dec = np.meshgrid(ra, dec)
            x, y = ra.ravel(), dec.ravel()
            img = healpy.get_interp_val(
                map_theta_phi, y * np.pi / 180.0, x * np.pi / 180.0
            )
        elif projection == "Robinson":
            ra = x = np.linspace(-180, 180, nside)
            dec = y = np.linspace(-90, 90, nside)
            ra, dec = np.meshgrid(ra, dec)

            img = healpy.get_interp_val(map_theta_phi, ra, dec, lonlat=True)
        else:
            return map_theta_phi

        img = img.reshape(nside, nside)

        return img

    def plot_map(
        self,
        freq=None,
        img=None,
        num=1,
        fig=None,
        ax=None,
        fig_kwargs={},
        include_cbar=True,
        projection="rectilinear",
        **kwargs,
    ):
        """
        Plot map in some coordinate system.

        .. note :: If the map is binary, e.g., a survey footprint, and the
            user supplies the keyword argument 'hatches', will overplot
            survey area is cross-hatched region (currently for 'rectilinear'
            projection only).

        .. note :: Can supply `img` directly. Will retrieve automatically
            if not supplied using `get_map` (see above).

        Parameters
        ----------
        freq : int, float
            Currently only applies if we're dealing with the GSM, in which case
            it should be the desired frequency in Hz.
        projection : str
            Desired projection of plot, options include 'rectlinear',
            'ecliptic', 'galactic', or 'equatorial' at the moment.

        Returns
        -------
        tuple
            A tuple containing the (figure object, axis object, projection
            object, and the image itself).
        """
        # Get image if it wasn't passed
        if img is None:
            if type(self.data) == str:
                if self.data == "gsm":
                    assert freq is not None, "Must supply `img` or `freq`!"

            img = self.get_map(freq=freq, projection=projection)

        # Setup plot window
        has_ax = True
        if ax is None:
            fig = plt.figure(num=num, **fig_kwargs)
            if projection.lower() in mollview_projections:
                # rect = 0, 1, 0, 1
                # ax = healpy.projaxes.MollweideAxes(fig=fig, rect=rect)
                has_ax = True
            else:
                has_ax = False

        # Setup kwargs for imshow or mollview call.
        if projection.lower() in mollview_projections:
            _kwargs = {}
            _proj = None
        elif projection.lower() == "robinson":
            try:
                import cartopy.crs as ccrs
            except ImportError:
                raise ImportError("Must have cartopy installed for Robinson projection")

            _proj = ccrs.Robinson(central_longitude=110)
            _kwargs = {
                "transform": ccrs.PlateCarree(),
                "extent": [-180, 180, -90, 90],
                "aspect": None,
                "origin": "lower",
            }
        elif projection == "rectilinear":
            _proj = None
            _kwargs = {"extent": [0, 24, -90, 90], "aspect": "auto"}

        kw = kwargs.copy()
        kw.update(_kwargs)

        if not has_ax:
            ax = fig.add_subplot(111, projection=_proj)

        ##
        # Actually plot
        if projection in mollview_projections:
            healpy.mollview(img, **kw)
            img = None
        elif projection == "rectilinear":
            if self.is_binary and "hatches" in kw:
                if type(kw["hatches"]) == str:
                    tmp = ["", kw["hatches"]]
                    kw["hatches"] = tmp

                del kw["aspect"]
                kw["origin"] = "image"
                ax.contourf(img, levels=1, **kw)
            else:
                ax.imshow(img, **kw)

            ax.set_xlim(24, 0)
            ax.set_ylim(-90, 90)
            ax.tick_params(direction="in", color="w", size=6, labelsize=22)
            ax.set_xlabel(r"Right Ascension [hours]", fontsize=24, labelpad=5)
            ax.set_ylabel(r"Declination [deg]", fontsize=24, labelpad=5)
        else:
            raise NotImplementedError("help")

        # May need transform for subsequent over-plotting
        proj = kwargs["transform"] if "transform" in kwargs else None

        return fig, ax, proj, img
