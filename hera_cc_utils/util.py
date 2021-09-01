""" Random stuff. """

import numpy as np
import healpy as hp

deg_per_hr = 15.

def field_to_healpix(dec_min, dec_max, ra_min=0, ra_max=24, nside=512):
    """
    Create a healpix map from a contiguous block in RA and DEC.


    """

    # Convert nside to number of pixels
    npix = hp.nside2npix(nside)

    # Create angles in radians
    theta, phi = hp.pix2ang(nside, np.arange(npix))

    # Convert to RA and DEC in hours and degrees
    ra = np.rad2deg(phi) / deg_per_hr
    dec = np.rad2deg(0.5 * np.pi - theta)

    # Setup a mask that tells us which pixels are in the HERA stripe.
    dec_ok = np.logical_and(dec >= dec_min, dec <= dec_max)
    ra_ok = np.logical_and(ra >= ra_min, ra <= ra_max)
    img = np.logical_and(dec_ok, ra_ok)

    return img

def get_overlap_area(img1, img2):
    """
    Given two healpix maps, find the overlapping area in square degrees.

    .. note :: Meant to be used with survey footprints, for which the maps are
        just arrays of ones and zeros, with ones indicating pixels covered
        by the survey.

    Parameters
    ----------
    img1, img2: np.ndarray
        Healpix maps (i.e., 1-D arrays).


    Returns
    -------
    Overlapping area in square degrees.
    """

    nside = hp.npix2nside(img1.size)
    nside2 = hp.npix2nside(img2.size)

    assert nside == nside2, "Maps must be at same resolution!"

    # Each pixel is this many square degrees
    pix = hp.nside2pixarea(nside, degrees=True)

    # Convert to boolean array
    overlap = (img1 > 0) * (img2 > 0)

    return overlap.sum() * pix

def get_map_area(img):
    """
    Given a healpix map, `img`, find the total area in square degrees.

    .. note :: Meant to be used with survey footprints, for which the 'map' is
        just an array of ones and zeros, with ones indicating pixels covered
        by the survey.

    Parameters
    ----------
    img : np.ndarray
        Healpix map.

    Returns
    -------
    Survey area in square degrees.

    """
    nside = hp.npix2nside(img.size)
    pix = hp.nside2pixarea(nside, degrees=True)

    mask = np.zeros_like(img)
    mask[img > 0] = 1

    return pix * img.sum()
