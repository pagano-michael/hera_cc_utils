"""Simple demo showing how to plot GSM."""

import healpy
import numpy as np
import hera_cc_utils as hera_cc

proj = "rectilinear"

gsm = hera_cc.Map()
rst = hera_cc.Map(data="roman")

kw = {"vmin": 1e-2} if proj == "rectilinear" else {"min": 1e-2}

fig1, ax1, proj1, img1 = gsm.plot_map(freq=150e6, projection=proj, num=1)

# Plot
fig2, ax2, proj2, img2 = rst.plot_map(projection=proj, num=2, cmap="binary", **kw)


# Just pick same nside as GSM
nside = 512
npix = healpy.nside2npix(nside)

# angles in radians
theta, phi = healpy.pix2ang(nside, np.arange(npix))

deg_per_hr = 15.0
ra = np.rad2deg(phi) / deg_per_hr
dec = np.rad2deg(0.5 * np.pi - theta)

hera_stripe = np.logical_and(dec >= -35, dec <= -25)

img = np.zeros(npix)
img[hera_stripe == 1] = 1

hera = hera_cc.Map(data=img)
fig3, ax3, proj3, img3 = hera.plot_map(projection=proj, num=3, cmap="binary_r", **kw)


img2[img2 < 1e-1] = 0

ax1.axhline(-35, color="w")
ax1.axhline(-25, color="w")

_kwargs = {"extent": [-6, 18, -90, 90], "aspect": "auto"}
ax1.imshow(img2, cmap="binary", alpha=0.4, **_kwargs)

lst_colors = ["limegreen", "gold", "cyan"]
lst_cuts = [(1.25, 2.7), (4.5, 6.5), (8.5, 10.75)]
for i, lc in enumerate(lst_cuts):
    ax1.fill_between(
        [lc[0], lc[1]],
        -34.5,
        -25.5,
        facecolor="none",
        edgecolor=lst_colors[i],
        alpha=1,
        zorder=5,
        lw=1.5,
    )
    ax1.text(lc[1] - 0.1, -32, "F{}".format(i + 1), fontsize=18, c="w")
