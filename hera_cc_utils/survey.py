""" Utilities for comparing k-space covered by different surveys.
    Thanks, Adrian! """

import numpy as np
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

class Survey:
    def __init__(self, survey_name, target_lines=None,
                 freq_range=None, lambda_range=None,
                 R=None, hemisphere=None, angular_res=None):
        """
        Initialize a survey
        """
        self.survey_name = survey_name

        if freq_range == None and lambda_range == None:
            raise ValueError('Must specify either wavelength or frequency range')
        else:
            if freq_range == None:
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
        self.R = R
        kpara_max = self.calc_kpara_max()
        kpara_min = self.calc_kpara_min()
        self.kpara_range = {}
        for line in self.target_lines:
            line_name = line.get_line_name()
            self.kpara_range[line_name] = \
            np.array([kpara_min[line_name], kpara_max[line_name]])
        kperp_max = self.calc_kperp_max()
        self.kperp_range = {}
        for line in self.target_lines:
            line_name = line.get_line_name()
            self.kperp_range[line_name] = \
            np.array([0.01, kperp_max[line_name]])

    def get_kperp_range(self):
        return self.kperp_range

    def get_hemisphere(self):
        return self.hemisphere

    def get_kpara_range(self):
        return self.kpara_range

    def calc_kperp_max(self):
        kperp_max = {}
        base_factors = np.pi / (self.angular_res * cosmo.h)
        for line in self.target_lines:
            min_z, max_z = self.redshift_ranges[line.get_line_name()]
            kperp_max[line.get_line_name()] = \
            base_factors / cosmo.comoving_distance(max_z).to_value()
        return kperp_max

    def calc_kpara_max(self):
        kpara_max = {}
        base_factors = 10**3 * self.R * np.pi * cosmo.H0.to_value() \
                        / (const.c.to_value() * cosmo.h)
        for line in self.target_lines:
            min_z, max_z = self.redshift_ranges[line.get_line_name()]
            kpara_max[line.get_line_name()] = \
            base_factors * cosmo.efunc(min_z) / (1 + min_z)
        return kpara_max

    def calc_kpara_min(self):
        kpara_min = {}
        for line in self.target_lines:
            min_z, max_z = self.redshift_ranges[line.get_line_name()]
            rpara_range = cosmo.comoving_distance(max_z) - cosmo.comoving_distance(min_z)
            rpara_range = rpara_range.to_value()
            kpara_min[line.get_line_name()] = 2. * np.pi / ( rpara_range * cosmo.h)
        return kpara_min

    def get_survey_name(self):
        return self.survey_name

    def get_target_lines(self):
        return [line.get_line_name() for line in self.target_lines]

    def get_freq_range(self):
        return self.freq_range

    def get_lambda_range(self):
        return self.lambda_range

    def calc_redshift_ranges(self):
        """
        Calculate a dictionary with the line name and the redshift range
        """
        redshift_ranges = {}
        for line in self.target_lines:
            min_redshift = np.max([line.freq_to_z(np.max(self.freq_range)), 0.])
            max_redshift = np.max([line.freq_to_z(np.min(self.freq_range)), 0.])
            redshift_ranges[line.get_line_name()] = np.array([min_redshift, max_redshift])
        return redshift_ranges

    def get_redshift_ranges(self,lines='All'):
        if lines == 'All':
            return self.redshift_ranges
        else:
            redshift_ranges = {}
            for line in lines:
                redshift_ranges[line.get_line_name()] = self.redshift_ranges[line.get_line_name()]
            return redshift_ranges
