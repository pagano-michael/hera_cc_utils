# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Simple demo showing how to plot GSM."""

import hera_cc_utils as hera_cc

gsm = hera_cc.Map()

fig1, ax1, proj1 = gsm.plot_map(freq=150e6, projection="Robinson")

fig2, ax2, proj2 = gsm.plot_map(
    freq=150e6, projection="rectilinear", num=2, vmin=1e2, vmax=1e3
)
