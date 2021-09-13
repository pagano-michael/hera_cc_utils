# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the MIT License

"""Utilities for dealing with common emission lines."""

from astropy import constants as const


class Line(object):
    """
    Define an object for handling emission lines.

    Parameters
    ----------
    line_name : str
        The name of the line.
    freq_rest : float
        The rest frame emission frequency, in MHz. Cannot be specified if
        `lambda_rest` is specified.
    lambda_rest : float
        The rest frame emission wavelength, in µm. Cannot be specified if
        `freq_rest` is specified.

    Raises
    ------
    ValueError
        This is raised if both `freq_rest` and `lambda_rest` are specified, or
        if neither is specified.
    """

    def __init__(self, line_name, freq_rest=None, lambda_rest=None):
        """Initialize with frequencies in MHz or wavelengths in microns."""
        if freq_rest is None and lambda_rest is None:
            raise ValueError("Must specify either rest wavelength or frequency")
        if freq_rest is not None and lambda_rest is not None:
            raise ValueError("Cannot specify both freq_rest and lambda_rest")

        self.line_name = line_name

        if freq_rest is None:
            self.lambda_rest = lambda_rest
            # Dividing m/s by microns gives frequencies in MHz
            self.freq_rest = const.c.to_value() / self.lambda_rest
        else:
            self.freq_rest = freq_rest
            # Dividing m/s by MHz gives wavelengths in microns
            self.lambda_rest = const.c.to_value() / self.freq_rest

    def get_line_name(self):
        """
        Get the name of the line.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The line name.
        """
        return self.line_name

    def get_rest_lambda(self):
        """
        Get the rest wavelength of the line.

        Parameters
        ----------
        None

        Returns
        -------
        float
            The rest wavelength of the line, in µm.
        """
        return self.lambda_rest

    def get_rest_freq(self):
        """
        Get the rest frequency of the line.

        Parameters
        ----------
        None

        Returns
        -------
        float
            The rest frequency of the line, in MHz.
        """
        return self.freq_rest

    def z_to_lambda(self, z):
        """
        Return the observed wavelength given the redshift.

        Parameters
        ----------
        z : float
            The redshift corresponding to the observation.

        Returns
        -------
        float
            The wavelength of the observed radiation, in µm.
        """
        return self.lambda_rest * (1.0 + z)

    def z_to_freq(self, z):
        """
        Return the observed frequency given the redshift.

        Parameters
        ----------
        z : float
            The redshift corresponding to the observation

        Returns
        -------
        float
            The frequency of the observed radiation, in MHz.
        """
        return self.freq_rest / (1.0 + z)

    def lambda_to_z(self, lambda_obs):
        """
        Return the redshift given the observed wavelength.

        Parameters
        ----------
        lambda_obs : float
            The observed wavelength of the radiation, in µm.

        Returns
        -------
        float
            The redshift corresponding to the observation.
        """
        return lambda_obs / self.lambda_rest - 1.0

    def freq_to_z(self, freq_obs):
        """
        Return the redshift given the observed frequency.

        Parameters
        ----------
        freq_obs : float
            The observed frequency of the radiation, in MHz.

        Returns
        -------
        float
            The redshift corresponding to the observation.
        """
        return self.freq_rest / freq_obs - 1.0
