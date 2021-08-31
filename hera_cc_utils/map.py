""" Utilities for dealing with maps. """

import os
import sys
import numpy as np
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle, Circle
from matplotlib.transforms import IdentityTransform
#from hera_cc_utils.util import degree_ticks, top_cbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from input import PATH

try:
    import cartopy.crs as ccrs
except ImportError:
    pass

try:
    from pygsm import GlobalSkyModel
except ImportError:
    pass

try:
    import healpy
except ImportError:
    pass

deg_per_hr = 15.

# PATH = os.environ.get('HERA_CC_UTILS')

map_options = 'gsm', 'roman', 'ngrst', 'euclid_south', 'euclid_fornax'

_coords_in = \
{
 'gsm': 'G',
 'roman': 'E',
 'ngrst': 'E',
 'euclid_south' : 'C',
 'euclid_fornax' : 'C'
}

_filenames = \
{
 'roman': 'Roman_GRISM_mask.fits',
 'euclid_south': 'Euclid_deep_south.fits',
 'euclid_fornax': 'Euclid_deep_Fornax.fits',
}

def coords_in(data):
    if type(data) == str:
        if data in _coords_in:
            return _coords_in[data]

    print("Assuming input map is in celestial coordinates.")
    return 'C'

class Map(object):
    """
    Utilities for dealing with maps, particularly (for now at least) all-sky
    maps like the global sky model (GSM).

    Parameters
    ----------
    data : str
        Data source for map. Current options include 'gsm'.
    freq_unit : str
        Assumed frequency units when generating maps. Options: Hz, MHz, GHz.
    kwargs : dict
        Any other optional keyword arguments accepted by the map being
        generated, e.g., that can be supplied to pygsm.GlobalSkyModel.
    """
    def __init__(self, data='gsm', freq_unit='Hz', coords_in=None, **kwargs):
        self.data = data
        self.freq_unit = freq_unit
        self.kwargs = kwargs

        if type(data) == str:
            assert self.data in map_options

    @property
    def _gsm(self):
        """ An instance of pygsm.GlobalSkyModel. """
        if not hasattr(self, '_gsm_'):
            self._gsm_ = GlobalSkyModel(freq_unit=self.freq_unit, **self.kwargs)
        return self._gsm_

    def get_map(self, freq=150e6, projection='ecliptic'):
        """
        Get a map suitable for plotting.

        .. note :: Currently this isn't general at all, just designed to work
            with GSM. Could be generalized to, e.g., plot Planck maps in the
            background.

        """

        coord_in = coords_in(self.data)

        if self.data in ['gsm']:
            raw_map = self._gsm.generate(freq)
            nside = 512
            rot = healpy.Rotator(coord=[coord_in, 'C'], deg=False, inv=True)
        elif self.data in ['roman', 'ngrst']:
            fn = _filenames['roman']
            raw_map = healpy.fitsfunc.read_map('{}/roman/{}'.format(PATH,
                fn))

            nside = healpy.npix2nside(raw_map.size)
            rot = healpy.Rotator(coord=[coord_in, 'C'], deg=False, inv=True)
        elif self.data in ['euclid_south','euclid_fornax']:
            fn = _filenames[self.data]
            raw_map = healpy.fitsfunc.read_map('{}/euclid/{}'.format(PATH,
                fn))
            nside = healpy.npix2nside(raw_map.size)
            rot = healpy.Rotator(coord=[coord_in, 'G'], deg=False, inv=True)
        elif type(self.data) == np.ndarray:
            raw_map = self.data
            nside = healpy.npix2nside(raw_map.size)
            rot = None
        else:
            # Could add other maps here for backdrops at some point.
            raise NotImplemented('help')

        ##
        # If we're here, need to do some transformations.

        theta, phi = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)))

        # Rotate into celestial coordinates
        if rot is not None:
            theta_cel, phi_cel = rot(theta, phi)
            map_theta_phi = healpy.get_interp_val(raw_map, theta_cel, phi_cel)
        else:
            theta_cel, phi_cel = theta, phi
            map_theta_phi = raw_map

        # setup RA Dec coords


        # interpolate healpix
        if projection == 'rectilinear':

            # Interpolate onto a rectangular mesh with size (nside, nside)
            ra = np.linspace(-90, 270, nside)
            dec = np.linspace(0, 180, nside)

            RA, DEC = np.meshgrid(ra, dec)
#
            X, Y = RA.ravel(), DEC.ravel()
            img = healpy.get_interp_val(map_theta_phi,
                Y * np.pi / 180., X * np.pi / 180.)
        elif projection == 'Robinson':
            ra = x = np.linspace(-180, 180, nside)
            dec = y = np.linspace(-90, 90, nside)
            RA, DEC = np.meshgrid(ra, dec)

            img = healpy.get_interp_val(map_theta_phi, RA, DEC, lonlat=True)
        else:
            return raw_map
            #_ra = np.rad2deg(phi) / deg_per_hr
            #_dec = np.rad2deg(0.5 * np.pi - theta)
            #ra = np.unique(_ra)
            #dec = np.unique(_dec)
#
            #RA, DEC = np.meshgrid(ra, dec)
#
            #img = healpy.get_interp_val(map_theta_phi, RA, DEC, lonlat=True)


        img = img.reshape(nside, nside)

        return img

    def plot_map(self, freq=None, img=None, num=1, fig=None, ax=None,
        fig_kwargs={}, include_cbar=True, projection='rectilinear', **kwargs):
        """
        Plot map in some coordinate system

        .. note :: Can supply `img` directly. Will retrieve automatically
            if not supplied using `get_map` (see above).

        Parameters
        -
        """

        coord_in = coords_in(self.data)

        # Get image if it wasn't passed
        if img is None:
            if self.data == 'gsm':
                assert freq is not None, "Must supply `img` or `freq`!"

            img = self.get_map(freq=freq, projection=projection)

        # Setup plot window
        has_ax = True
        if ax is None:
            fig = plt.figure(num=num, **fig_kwargs)
            if projection.lower() in ['ecliptic', 'galactic']:
                rect = 0, 1, 0, 1
                #ax = healpy.projaxes.MollweideAxes(fig=fig, rect=rect)
                has_ax = True
            else:
                has_ax = False

        # Setup kwargs for imshow or mollview call.
        if projection.lower() in ['ecliptic', 'galactic']:
            _kwargs = {}
            _proj = None
        elif projection.lower() == 'robinson':
            _proj = ccrs.Robinson(central_longitude=110)
            _kwargs = {'transform': ccrs.PlateCarree(),
                'extent': [-180, 180, -90, 90], 'aspect': None,
                'origin':'lower'}
        elif projection == 'rectilinear':
            _proj = None
            _kwargs = {'extent': [-6, 18, -90, 90], 'aspect': 'auto'}

        kw = kwargs.copy()
        kw.update(_kwargs)

        if not has_ax:
            ax = fig.add_subplot(111, projection=_proj)

        ##
        # Actually plot
        if projection == 'ecliptic':
            healpy.mollview(img, coord=[coord_in, 'E'], **kw)
            img = None
        elif projection == 'galactic':
            healpy.mollview(img, coord=[coord_in, 'G'], **kw)
            img = None
        elif projection == 'rectilinear':
            cax = ax.imshow(img, **kw)

            ax.set_xlim(18, -6)
            ax.set_ylim(-90, 90)
            ax.tick_params(direction='in', color='w', size=6, labelsize=22)
            ax.set_xlabel(r'Right Ascension [hours]', fontsize=24, labelpad=5)
            ax.set_ylabel(r'Declination [deg]', fontsize=24, labelpad=5)
            fig.tight_layout()
        else:
            raise NotImplemented('help')

        # May need transform for subsequent over-plotting
        proj = kwargs['transform'] if 'transform' in kwargs else None

        return fig, ax, proj, img
