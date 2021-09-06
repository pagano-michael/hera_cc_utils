""" Utilities for comparing k-space covered by different surveys.
    Thanks, Adrian! """

import numpy as np
from .line import Line
import matplotlib.pyplot as plt
from .util import line_registry
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

survey_registry = \
{
 'hera': {'target_lines': ['HI'], 'R': None, 'hemisphere': 'S',
    'freq_range': [50.,200.], 'angular_res': None},

 'spherex': {'target_lines': ['Ha', 'Hb', 'Lya', 'OII', 'OIII'],
    'R': 150., 'hemisphere': 'B',
    'lambda_range': [0.75,5.0], 'angular_res':3.0058*10**-5},

 'fyst': {'target_lines': ['CII', 'CO21', 'CO32', 'CO43', 'CO54'],
   'R': 100., 'hemisphere': 'S',
   'freq_range': [220000., 405000.], 'angular_res':0.000193925},

 'tim': {'target_lines': ['CII', 'CO32', 'CO43', 'CO54'],
   'R': 250., 'hemisphere': 'S',
   'freq_range': [240000., 420000.], 'angular_res':5.8178*10**-5},

 'concerto': {'target_lines': ['CII', 'CO21', 'CO32', 'CO43', 'CO54'],
   'R': 200., 'hemisphere': 'S',
   'freq_range': [200000., 360000.], 'angular_res':5.8178*10**-5},
}

survey_registry['ccatp'] = survey_registry['fyst']

class Survey(object):
    def __init__(self, survey_name, target_lines=None,
                 freq_range=None, lambda_range=None,
                 R=None, hemisphere=None, angular_res=None):
        """
        Initialize a survey
        """
        self.survey_name = survey_name

        if survey_name in survey_registry.keys():

            reg = survey_registry[survey_name]

            # More elegant way to do this but I'm lazy
            if 'freq_range' in reg:
                freq_range = reg['freq_range']
            if 'lambda_range' in reg:
                lambda_range = reg['lambda_range']
            if 'angular_res' in reg:
                angular_res = reg['angular_res']
            if 'target_lines' in reg:
                target_lines = [Line(line, **line_registry[line]) \
                    for line in reg['target_lines']]
            if 'R' in reg:
                R = reg['R']

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

    #def get_line_wave(self):
    #
    #    return [line.get_line_name() for line in self.target_lines]

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

    def plot_coverage_k(self, ax=None, fig=1, **kwargs):
        """
        Plot a rectangle in (k_perp, k_para) space for each line targeted
        by the survey.
        """

        had_ax = True
        if ax is None:
            fig, ax = plt.subplots(1, 1, num=fig)
            had_ax = False

        kpara = self.get_kpara_range()
        kperp = self.get_kperp_range()

        for line in self.get_target_lines():
            ax.fill_between(kperp[line], *kpara[line], **kwargs)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-3, 1e2)
        ax.set_ylim(1e-3, 1e2)

        if not had_ax:
            ax.set_xlabel(r'$k_\perp$ [$h$Mpc$^{-1}$]')
            ax.set_ylabel(r'$k_\parallel$ [$h$Mpc$^{-1}$]')

        return ax

    def plot_coverage_z(self, ax=None, fig=1, use_labels=True, start=0,
        **kwargs):
        """
        Plot a horizontal line for each emission line targeted by survey vs.
        redshift on the x-axis.
        """

        had_ax = True
        if ax is None:
            fig, ax = plt.subplots(1, 1, num=fig)
            had_ax = False

        line_zs = self.get_redshift_ranges()
        line_names = self.get_target_lines()

        for i, line in enumerate(line_names):
            label = '{} {}'.format(self.survey_name, line) if use_labels \
                else None
            ax.plot(line_zs[line], [start+i]*2, label=label, **kwargs)

        ax.set_xlim(0, 15)
        ax.set_yticklabels([])
        ax.legend()

        if not had_ax:
            ax.set_xlabel(r'$z$')

        return ax
