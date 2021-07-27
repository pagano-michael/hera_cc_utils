""" Utilities for deailng with maps. """

import sys
import numpy as np
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle, Circle
from matplotlib.transforms import IdentityTransform
#from hera_cc_utils.util import degree_ticks, top_cbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    import cartopy.crs as ccrs
except ImportError:
    pass

try:
    from pygsm import GlobalSkyModel
except IndexError:
    pass

try:
    import healpy
except ImportError:
    pass

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
    def __init__(self, data='gsm', freq_unit='Hz', **kwargs):
        self.data = data
        self.freq_unit = freq_unit
        self.kwargs = kwargs

        assert self.data in ['gsm']

    @property
    def _gsm(self):
        """ An instance of pygsm.GlobalSkyModel. """
        if not hasattr(self, '_gsm_'):
            self._gsm_ = GlobalSkyModel(freq_unit=self.freq_unit, **self.kwargs)
        return self._gsm_

    def get_map(self, freq=150e6, nside=512, projection='rectilinear'):
        """
        Get a map suitable for plotting.

        .. note :: Currently this isn't general at all, just designed to work
            with GSM. Could be generalized to, e.g., plot Planck maps in the
            background.

        """
        if self.data in ['gsm']:
            raw_map = self._gsm.generate(freq)
        else:
            # Could add other maps here for backdrops at some point.
            raise NotImplemented('help')

        theta, phi = healpy.pix2ang(nside, np.arange(healpy.nside2npix(nside)))

        # rotate into celestial
        rot = healpy.Rotator(coord=['G', 'C'], deg=False, inv=True)
        theta_cel, phi_cel = rot(theta, phi)
        map_theta_phi = healpy.get_interp_val(raw_map, theta_cel, phi_cel)

        # setup RA Dec coords


        # interpolate healpix
        if projection == 'rectilinear':
            nside = 1024
            ra = np.linspace(-90, 270, nside)
            dec = np.linspace(0, 180, nside)
            RA, DEC = np.meshgrid(ra, dec)

            X, Y = RA.ravel(), DEC.ravel()
            img = healpy.get_interp_val(map_theta_phi,
                Y * np.pi / 180., X * np.pi / 180.)
        else:
            ra = x = np.linspace(-180, 180, nside)
            dec = y = np.linspace(-90, 90, nside)
            RA, DEC = np.meshgrid(ra, dec)

            img = healpy.get_interp_val(map_theta_phi, RA, DEC, lonlat=True)


        img = img.reshape(nside, nside)

        return img

    def plot_map(self, freq=None, nside=512, img=None, num=1, fig=None, ax=None,
        fig_kwargs={}, include_cbar=True, projection='rectilinear', **kwargs):
        """
        Plot map in some coordinate system

        .. note :: Can supply `img` directly. Will retrieve automatically
            if not supplied using `get_map` (see above).

        Parameters
        -
        """

        # Get image if it wasn't passed
        if img is None:
            assert freq is not None, "Must supply `img` or `freq`!"
            img = self.get_map(freq, nside, projection)

        # Setup plot window
        has_ax = True
        if ax is None:
            fig = plt.figure(num=num, **fig_kwargs)
            has_ax = False

        if projection == 'Robinson':
            _proj = ccrs.Robinson(central_longitude=110)
            _kwargs = {'transform': ccrs.PlateCarree(),
                'extent': [-180, 180, -90, 90], 'aspect': None,
                'origin':'lower'}
        elif projection == 'rectilinear':
            _proj = None
            _kwargs = {'extent': [-6, 18, -90, 90], 'aspect': 'auto',
                'origin':'lower'}

        kw = kwargs.copy()
        kw.update(_kwargs)

        if ('norm' not in kw) and ('vmin' not in kw) and ('vmax' not in kw):
            kw['norm'] = LogNorm()

        if not has_ax:
            ax = fig.add_subplot(111, projection=_proj)

        # Actually plot
        cax = ax.imshow(img, **kw)

        # Make nice axes
        if projection == 'rectilinear':
            ax.set_xlim(18, -6)
            ax.set_ylim(-90, 90)
            ax.tick_params(direction='in', color='w', size=6, labelsize=22)
            ax.set_xlabel(r'Right Ascension [hours]', fontsize=24, labelpad=5)
            ax.set_ylabel(r'Declination [deg]', fontsize=24, labelpad=5)


        # May need transform for subsequent over-plotting
        proj = kwargs['transform'] if 'transform' in kwargs else None

        return fig, ax, proj
