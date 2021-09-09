""" Utilities for dealing with common emission lines. Thanks, Adrian! """

from astropy import constants as const


class Line(object):
    def __init__(self, line_name, freq_rest=None, lambda_rest=None):
        """
        Initialize with frequencies in MHz or wavelengths in microns.
        """
        if freq_rest == None and lambda_rest == None:
            raise ValueError("Must specify either rest wavelength or frequency")
        else:
            self.line_name = line_name

            if freq_rest == None:
                self.lambda_rest = lambda_rest
                # Dividing m/s by microns gives frequencies in MHz
                self.freq_rest = const.c.to_value() / self.lambda_rest
            else:
                self.freq_rest = freq_rest
                # Dividing m/s by MHz gives wavelengths in microns
                self.lambda_rest = const.c.to_value() / self.freq_rest

    def get_line_name(self):
        return self.line_name

    def get_rest_lambda(self):
        """
        Return the rest wavelength
        """
        return self.lambda_rest

    def get_rest_freq(self):
        """
        Return the rest frequency
        """
        return self.freq_rest

    def z_to_lambda(self, z):
        """
        Return the observed wavelength given the redshift
        """
        return self.lambda_rest * (1.0 + z)

    def z_to_freq(self, z):
        """
        Return the observed frequency given the redshift
        """
        return self.freq_rest / (1.0 + z)

    def lambda_to_z(self, lambda_obs):
        """
        Return the redshift given the observed wavelength
        """
        return lambda_obs / self.lambda_rest - 1.0

    def freq_to_z(self, freq_obs):
        """
        Return the redshift given the observed frequency
        """
        return self.freq_rest / freq_obs - 1.0
