# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Utilities for comparing k-space covered by different surveys."""

import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from matplotlib.patches import Rectangle
from astropy.cosmology import FlatLambdaCDM

from .line import Line
from .util import line_registry

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

survey_registry = {
    "hera": {
        "target_lines": ["HI"],
        "rval": 1533.0,
        "hemisphere": "S",
        "freq_range": [50.0, 200.0],
        "angular_res": 0.00727,
    },
    "spherex": {
        "target_lines": ["Ha", "Hb", "Lya", "OII", "OIII"],
        "rval": 150.0,
        "hemisphere": "B",
        "lambda_range": [0.75, 5.0],
        "angular_res": 3.0058 * 1e-5,
    },
    "fyst": {
        "target_lines": ["CII", "CO21", "CO32", "CO43", "CO54"],
        "rval": 100.0,
        "hemisphere": "S",
        "freq_range": [220000.0, 405000.0],
        "angular_res": 0.000193925,
    },
    "tim": {
        "target_lines": ["CII", "CO32", "CO43", "CO54"],
        "rval": 250.0,
        "hemisphere": "S",
        "freq_range": [240000.0, 420000.0],
        "angular_res": 5.8178 * 1e-5,
    },
    "concerto": {
        "target_lines": ["CII", "CO21", "CO32", "CO43", "CO54"],
        "rval": 200.0,
        "hemisphere": "S",
        "freq_range": [200000.0, 360000.0],
        "angular_res": 5.8178 * 1e-5,
    },
    "roman": {
        "target_lines": ["Lya"],
        "rval": 461, # Really 461 * wave, with `wave` in microns.
        "hemisphere": "S",
        "lambda_range": [1, 1.93],
        "angular_res": 5.333e-07, # 0.11 arcsec/pixel
    },
    "euclid": {
        "target_lines": ["Lya"],
        "rval": 250,
        "hemisphere": "S",
        "lambda_range": [0.92, 2.],
        "angular_res": 4.848e-07, # 0.1 arcsec/pixel
    }
}

survey_registry["ccatp"] = survey_registry["fyst"]


class Survey(object):
    """
    Define an object for handling surveys.

    Parameters
    ----------
    survey_name : str
        The name of the survey. List of known surveys is: hera, spherex, fyst,
        tim, concerto.
    target_lines : list of str
        The lines to include in the survey. Ignored if `survey_name` matches one
        of the known surveys.
    freq_range : 2-element list of float
        The range of frequencies covered by the survey, in MHz. Ignored if
        `survey_name` matches one of the known surveys. Should not be specified
        if `lambda_range` is specified.
    lambda_range : 2-element list of float
        The range of wavelengths covered by the survey, in µm. Ignored if
        `survey_name` matches one of the known surveys. Should not be specified
        if `freq_range` is specified.
    rval : float
        The resolving power "R" of the experiment. Ignored if `survey_name`
        matches one of the known surveys.
    hemisphere : str
        The hemisphere the survey measures. Should be one of: "N" (north), "S"
        (south), or "B" (both).
    angular_res : float
        The angular resolution of the experiment, in radians. Ignored if
        `survey_name` matches one of the known surveys.

    Raises
    ------
    ValueError
        This is raised if `survey_name` is not one of the known surveys, and
        either both `freq_range` and `lambda_range` are specified or neither is
        specified.

    """

    def __init__(
        self,
        survey_name,
        target_lines=None,
        freq_range=None,
        lambda_range=None,
        rval=None,
        hemisphere=None,
        angular_res=None,
    ):
        """Initialize a survey."""
        self.survey_name = survey_name

        if survey_name in survey_registry.keys():
            reg = survey_registry[survey_name]

            # More elegant way to do this but I'm lazy
            if "freq_range" in reg:
                freq_range = reg["freq_range"]
            if "lambda_range" in reg:
                lambda_range = reg["lambda_range"]
            if "angular_res" in reg:
                angular_res = reg["angular_res"]
            if "target_lines" in reg:
                target_lines = [
                    Line(line, **line_registry[line]) for line in reg["target_lines"]
                ]
            if "rval" in reg:
                rval = reg["rval"]

        if freq_range is None and lambda_range is None:
            raise ValueError("Must specify either wavelength or frequency range")
        if freq_range is not None and lambda_range is not None:
            raise ValueError("Cannot specify freq_range and lambda_range")
        else:
            if freq_range is None:
                self.lambda_range = np.array(lambda_range)
                # Dividing m/s by microns gives frequencies in MHz
                self.freq_range = const.c.to_value() / self.lambda_range
            else:
                self.freq_range = np.array(freq_range)
                # Dividing m/s by MHz gives wavelengths in microns
                self.lambda_range = const.c.to_value() / self.freq_range

        self.angular_res = angular_res
        self.hemisphere = hemisphere
        self.target_lines = target_lines
        self.redshift_ranges = self.calc_redshift_ranges()
        self.rval = rval
        kpara_max = self.calc_kpara_max()
        kpara_min = self.calc_kpara_min()
        self.kpara_range = {}
        for line in self.target_lines:
            line_name = line.get_line_name()
            self.kpara_range[line_name] = np.array(
                [kpara_min[line_name], kpara_max[line_name]]
            )
        kperp_max = self.calc_kperp_max()
        self.kperp_range = {}
        for line in self.target_lines:
            line_name = line.get_line_name()
            self.kperp_range[line_name] = np.array([0.01, kperp_max[line_name]])

    def get_kperp_range(self):
        """
        Get the k_perpendicular range covered by a survey.

        Parameters
        ----------
        None

        Returns
        -------
        dict of 2-element list of float
            The minimum and maximum k_perpendicular values, in units of
            1/Mpc. The keys of the dictionary correspond to target_lines for the
            survey and the values are the lists of ranges.
        """
        return self.kperp_range

    def get_hemisphere(self):
        """
        Get the hemisphere covered by a survey.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The hemisphere covered by the survey.

        """
        return self.hemisphere

    def get_kpara_range(self):
        """
        Get the k_parallel range covered by a survey.

        Parameters
        ----------
        None

        Returns
        -------
        dict of 2-element list of float
            The minimum and maximum k_parallel values, in units of 1/Mpc. The
            keys of the dictionary correspond to target_lines for the survey and
            the values are the lists of ranges.
        """
        return self.kpara_range

    def calc_kperp_max(self):
        """
        Calculate the maximum k_perpendicular modes for a survey.

        Parameters
        ----------
        None

        Returns
        -------
        dict of float
            The minimum k_perpendicular value for a survey, in units of
            1/Mpc. The keys correspond to the target_lines of the survey.
        """
        kperp_max = {}
        base_factors = np.pi / (self.angular_res * cosmo.h)
        for line in self.target_lines:
            min_z, max_z = self.redshift_ranges[line.get_line_name()]
            kperp_max[line.get_line_name()] = (
                base_factors / cosmo.comoving_distance(max_z).to_value()
            )
        return kperp_max

    def calc_kpara_max(self):
        """
        Calculate the maximum k_parallel modes for a survey.

        Parameters
        ----------
        None

        Returns
        -------
        dict of float
            The maximum k_parallel value for a survey, in units of 1/Mpc. The
            keys correspond to the target_lines of the survey.
        """
        kpara_max = {}
        base_factors = (
            10 ** 3
            * self.rval
            * np.pi
            * cosmo.H0.to_value()
            / (const.c.to_value() * cosmo.h)
        )
        for line in self.target_lines:
            min_z, max_z = self.redshift_ranges[line.get_line_name()]
            kpara_max[line.get_line_name()] = (
                base_factors * cosmo.efunc(min_z) / (1 + min_z)
            )
        return kpara_max

    def calc_kpara_min(self):
        """
        Calculate the minimum k_parallel modes for a survey.

        Parameters
        ----------
        None

        Returns
        -------
        dict of float
            The minimum k_parallel values for a survey, in units of 1/Mpc. The
            keys correspond to the target_lines of the survey.
        """
        kpara_min = {}
        for line in self.target_lines:
            min_z, max_z = self.redshift_ranges[line.get_line_name()]
            rpara_range = cosmo.comoving_distance(max_z) - cosmo.comoving_distance(
                min_z
            )
            rpara_range = rpara_range.to_value()
            kpara_min[line.get_line_name()] = 2.0 * np.pi / (rpara_range * cosmo.h)
        return kpara_min

    def get_survey_name(self):
        """
        Get the name of a survey.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The survey name.
        """
        return self.survey_name

    def get_target_lines(self):
        """
        Get the targeted transition lines for a survey.

        Parameters
        ----------
        None

        Returns
        -------
        list of str
            The names of the lines.
        """
        return [line.get_line_name() for line in self.target_lines]

    def get_freq_range(self):
        """
        Get the frequency range for a survey.

        Parameters
        ----------
        None

        Returns
        -------
        2-element list of float
            The minimum and maximum frequencies of the survey, in MHz.
        """
        return self.freq_range

    def get_lambda_range(self):
        """
        Get the wavelength range for a survey.

        Parameters
        ----------
        None

        Returns
        -------
        2-element list of float
            The minimum and maximum wavelengths of the survye, in µm.
        """
        return self.lambda_range

    def calc_redshift_ranges(self):
        """
        Calculate a dictionary with the line name and the redshift range.

        Parameters
        ----------
        None

        Returns
        -------
        dict of 2-element 1d-arrays of float
            The minimum and maximum redshift values for a survey. Dictionary
            keys correspond to target_lines, and values are the redshift range.
        """
        redshift_ranges = {}
        for line in self.target_lines:
            min_redshift = np.max([line.freq_to_z(np.max(self.freq_range)), 0.0])
            max_redshift = np.max([line.freq_to_z(np.min(self.freq_range)), 0.0])
            redshift_ranges[line.get_line_name()] = np.array(
                [min_redshift, max_redshift]
            )
        return redshift_ranges

    def get_redshift_ranges(self, lines=None):
        """
        Get the redshift range for a survey for a particular line.

        Parameters
        ----------
        lines : list of str, optional
            The list of lines to get redshifts for. If `None`, then all redshift
            ranges are returned.

        Returns
        -------
        dict of 2-element 1d-arrays of float
            The minimum and maximum redshift values for a survey. Dictionary
            keys correspond to target_lines, and values are the redshift range.
        """
        if lines is None:
            return self.redshift_ranges
        else:
            redshift_ranges = {}
            for line in lines:
                redshift_ranges[line.get_line_name()] = self.redshift_ranges[
                    line.get_line_name()
                ]
            return redshift_ranges

    def plot_coverage_k(self, ax=None, fig=1, **kwargs):
        """
        Plot in (k_perp, k_para) space for each line targeted by the survey.

        Parameters
        ----------
        ax : matplotlib axis object, optional
            The axis to add the plot to. If None, create a new axis.
        fig : int, optional
            The figure number to attach to.
        kwargs : dict
            The options to pass through to `matplotlib.pyplot.fill_between`.

        Returns
        -------
        ax : matplotlib axis object
            The axis the figure was plotted to.
        """
        had_ax = True
        if ax is None:
            fig, ax = plt.subplots(1, 1, num=fig)
            had_ax = False

        kpara = self.get_kpara_range()
        kperp = self.get_kperp_range()

        for line in self.get_target_lines():
            rect = Rectangle((kperp[line][0], kpara[line][0]),
                np.diff(kperp[line])[0], np.diff(kpara[line])[0],
                **kwargs)
            ax.add_patch(rect)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-3, 1e2)
        ax.set_ylim(1e-3, 1e2)

        if not had_ax:
            ax.set_xlabel(r"$k_\perp$ [$h$Mpc$^{-1}$]")
            ax.set_ylabel(r"$k_\parallel$ [$h$Mpc$^{-1}$]")

        return ax

    def plot_coverage_z(self, ax=None, fig=1, use_labels=True, start=0, **kwargs):
        """
        Plot the redshift coverage for each line targeted by a survey.

        This method plots a horizontal line for each emission line targeted by
        survey vs.  redshift on the x-axis.

        Parameters
        ----------
        ax : matplotlib axis object, optional
            The axis to add the plot to. If None, create a new axis.
        fig : int, optional
            The figure number to attach to.
        use_labels : bool, optional
            Whether to add labels to the plot.
        start : int, optional
            Which line to start plotting with.
        kwargs : dict
            The keyword arguments passed to `matplotlib.pyplot.plot`.

        Returns
        -------
        ax : matplotlib axis object
            The axis the figure was plotted to.
        """
        had_ax = True
        if ax is None:
            fig, ax = plt.subplots(1, 1, num=fig)
            had_ax = False

        line_zs = self.get_redshift_ranges()
        line_names = self.get_target_lines()

        for i, line in enumerate(line_names):
            label = "{} {}".format(self.survey_name, line) if use_labels else None
            ax.plot(line_zs[line], [start + i] * 2, label=label, **kwargs)

        ax.set_xlim(0, 15)
        ax.set_yticklabels([])
        ax.legend()

        if not had_ax:
            ax.set_xlabel(r"$z$")

        return ax
