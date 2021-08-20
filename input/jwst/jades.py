"""
Information about JADES.

https://www.cosmos.esa.int/web/jwst-nirspec-gto/jades

So far, just program 1180:

https://www.stsci.edu/jwst/phase2-public/1180.pdf

particularly the deep fields.

"""

# Observation number in
jades_obs = [7,8,9,10,11,12,13,14,15,16,17,18]

# V3PA is the position angle (PA) of the V3 reference axis eastward relative to
# north when projected onto the sky.
v3pa = [295] * 6 + [299] * 6

# pre-imaging program GTO 1180
jades = \
    [('03h32m42.1560s', '-27d47m57.10s'),
     ('03h32m42.6680s', '-27d47m59.64s'),
     ('03h32m44.1300s', '-27d47m0.91s'),
     ('03h32m44.6420s', '-27d47m3.45s'),
     ('03h32m33.7260s', '-27d47m3.73s'),
     ('03h32m44.2380s', '-27d47m3.73s'), # end of deep pointings 1 & 2
     ('03h32m35.9910s', '-27d46m6.96s'),
     ('03h32m36.5030s', '-27d46m9.50s'),
     ('03h32m47.5440s', '-27d48m5.45s'),
     ('03h32m44.5410s', '-27d48m52.95s'),
     ('03h32m40.0960s', '-27d46m42.53s'),
     ('03h32m37.0930s', '-27d47m30.03s'), # end of deep pointings 3 & 4
     ]

# Position of NIRCAM on focal plane
# V2_Ref	V3_Ref	V3_IdlYAngle	V2_1	V2_2	V2_3	V2_4	V3_1	V3_2	V3_3	V3_4
nircam = -0.32, -492.59, -0.03, 153.16, -153.74, -152.07, 151.38, -559.27, \
    -557.14, -426.00, -427.95
