""" Random stuff. """

import numpy as np
import healpy as hp

deg_per_hr = 15.

def field_to_healpix(dec_min, dec_max, nside=512):

    # Convert nside to number of pixels
    npix = hp.nside2npix(nside)

    # Create angles in radians
    theta, phi = hp.pix2ang(nside, np.arange(npix))

    # Convert to RA and DEC in hours and degrees
    ra = np.rad2deg(phi) / deg_per_hr
    dec = np.rad2deg(0.5 * np.pi - theta)

    # Setup a mask that tells us which pixels are in the HERA stripe.
    img = np.logical_and(dec >= dec_min, dec <= dec_max)

    return img
