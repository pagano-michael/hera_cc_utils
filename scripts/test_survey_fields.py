""" Simple demo showing how to plot GSM. """

import numpy as np
import hera_cc_utils as hera_cc
from astropy import units as u
from astropy.coordinates import SkyCoord

gsm = hera_cc.Map()

fig1, ax1, proj1 = gsm.plot_map(freq=150e6, projection='Robinson')

fig2, ax2, proj2 = gsm.plot_map(freq=150e6, projection='rectilinear',
    num=2, vmin=1e2, vmax=1e3)

ax2.axhline(-30 - 5, color='w', ls=':', lw=1.5, zorder=1)
ax2.axhline(-30 + 5, color='w', ls=':', lw=1.5, zorder=1)

lst_cuts = [(1.25, 2.7), (4.5, 6.5), (8.5, 10.75)]
for i, lst_cut in enumerate(lst_cuts):
    ax2.fill_betweenx([-35, -25], *lst_cut, color='none', edgecolor='w')
    ax2.annotate('field {}'.format(i+1), (np.mean(lst_cut), -24), color='w',
        ha='center', va='bottom')

ax2.set_xlim(1, 11)
ax2.set_ylim(-40, -20)

# Draw GOODS-S field
goods_S = SkyCoord('03h32m28s', '-27d48m30s', frame='icrs')

ra, dec = goods_S.ra.hour, goods_S.dec.degree

ax2.scatter(ra, dec, marker='x', color='w')

#patch = Rectangle(xy=(ra,dec), width=0.5, height=0.5,
#    facecolor='c', alpha=0.5)
#ax.add_patch(patch)

euclid_dff = SkyCoord('03h31m43.6s', '-28d05m18.6s', frame='icrs')
euclid_sqdeg = 20
euclid_R = np.sqrt(10. / np.pi)

ra, dec = euclid_dff.ra.hour, euclid_dff.dec.degree
ax2.scatter(ra, dec, color='w', marker='+')
