#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018-05-14
# Last Modified: 2019-02-17
from __future__ import absolute_import, division, print_function
__author__     = ['Victor Calderon']
__copyright__  = ["Copyright 2018 Victor Calderon"]
__email__      = ['victor.calderon@vanderbilt.edu']
__maintainer__ = ['Victor Calderon']
__all__        = [  "Behroozi_relation",
                    "Moster2010_relation",
                    "Moster_example"]

## Import modules
import numpy as np
from   cosmo_utils.utils             import file_utils as fd
from   cosmo_utils.custom_exceptions import LSSUtils_Error

## Functions

## Retrieves default values for Behroozi et al. (2013)
def _retrieve_Behroozi_default_dict():
    """
    Dictionary of default values of all model parameters set to the
    column 2 values in Table 2 of Behroozi et al. (2013)

    Returns
    --------
    d : `dict`
        Dictionary containing default parameters for the Stellar-Halo
        Mass relation of Behroozi et al. (2013)

    Notes
    ----------
    All calculations are done internally ising the same h=0.7 units as
    in Behroozi ete al. (2010), ['arXiv:1001.0015'] so the parameter values
    here are the same as in Table 2, even though the `mean_log_halo_mass`
    and `mean_stellar_mass` methods use accept and return arguments in
    h=1 units.
    """
    ## Main dictionary
    d = ({'smhm_m0_0': 10.72,
        'smhm_m0_a': 0.59,
        'smhm_m1_0': 12.35,
        'smhm_m1_a': 0.3,
        'smhm_beta_0': 0.43,
        'smhm_beta_a': 0.18,
        'smhm_delta_0': 0.56,
        'smhm_delta_a': 0.18,
        'smhm_gamma_0': 1.54,
        'smhm_gamma_a': 2.52})

    return d

## Behroozi SHMR function
def Behroozi_relation(log_mstar, z=0.):
    """
    Returns the halo mass of a central galaxy as a function of its stellar
    mass.

    Parameters
    -----------
    log_mstar : `float` ,`np.ndarray`, or array-like
        Value or array of values of base-10 logarithm of stellar mass
        in h=1 solar mass units.

    z : int, float, `np.ndarray` or array-like
        Redshift of the halo hosting the galaxy. If passing an array,
        it must be of the same length as the input `log_mstar`.

    Returns
    -----------
    log_halo_mass : float or `np.ndarray`
        Array or float containing 10-base logarithm of halo mass in ``h=1``
        solar mass units.

    Notes
    ----------
    The parameter values in Behroozi+10 were fit to data assuming ``h=0.7``.
    Thus, we will transform our input stellar mass to ``h=0.7`` units,
    evaluate using the Behroozi parameters, and then transform back to
    ``h=1`` units before returning the result.
    """
    file_msg = fd.Program_Msg(__file__)
    little_h = 0.7
    ## Checking input parameters
    # `log_mstar`
    mstar_valid_types = (int, float, np.ndarray, list)
    if not (isinstance(log_mstar, mstar_valid_types)):
        msg = '{0} `log_mstar` ({1}) is not a valid type!'.format(
            file_msg, type(log_mstar))
        raise TypeError(msg)
    ##
    ## Behroozi dictionary
    param_dict = _retrieve_Behroozi_default_dict()
    ## COnverting stellar mass from ``h=1`` units to ``h=0.7`` units.
    mstar = (10.**log_mstar) / (little_h**2)
    ## Scale factor
    a = 1./(1. + z)
    ##
    ## Behroozi function
    logm0 = param_dict['smhm_m0_0'] + param_dict['smhm_m0_a']*(a - 1.)
    m0    = 10.**logm0
    logm1 = param_dict['smhm_m1_0'   ] + param_dict['smhm_m1_a'   ]*(a - 1)
    beta  = param_dict['smhm_beta_0' ] + param_dict['smhm_beta_a' ]*(a - 1)
    delta = param_dict['smhm_delta_0'] + param_dict['smhm_delta_a']*(a - 1)
    gamma = param_dict['smhm_gamma_0'] + param_dict['smhm_gamma_a']*(a - 1)
    #
    stellar_mass_by_m0 = mstar/m0
    term3_numerator    = (stellar_mass_by_m0)**delta
    term3_denominator  = 1. + (stellar_mass_by_m0)**(-gamma)
    #
    log_halo_mass = logm1 + beta*np.log10(stellar_mass_by_m0)
    log_halo_mass += (term3_numerator/term3_denominator) - 0.5
    #
    # Convert back from ``h=0.7`` to ``h=1`` units
    return np.log10((10.**log_halo_mass)*(little_h))

## Retrieves default values for Moster et al. (2010)
def _retrieve_Moster_default_dict():
    """
    Dictionary of default values of all model parameters set to the
    column 2 values in Table 2 of Moster et al. (2010)

    Returns
    --------
    d : `dict`
        Dictionary containing default parameters for the Stellar-Halo
        Mass relation of Moster et al. (2010)

    Notes
    ----------
    All calculations are done internally using the same h=0.7 units as
    in Moster et al. (2010), ['arXiv:0903.4682'] so the parameter values
    here are the same as in Table 2, even though the `mean_log_halo_mass`
    and `mean_stellar_mass` methods use accept and return arguments in
    h=1 units.
    """
    ## Main dictionary
    d = ({'smhm_m1':11.899,
        'smhm_m_M_0':0.02817,
        'smhm_beta':1.068,
        'smhm_gamma':0.611})

    return d

## Moster SHMR function
def Moster2010_relation(log_halo_mass, return_h0=False):
    """
    Returns the halo mass of a central galaxy as a function of its stellar
    mass.

    Parameters
    -----------
    log_halo_mass : `float` ,`np.ndarray`, or array-like
        Value or array of values of base-10 logarithm of halo mass
        in ``h = 1`` solar mass units.

    return_h0 : `bool`, optional
        If `True`, the function returns both stellar mass and halo mass
        arrays in units of ``h = 1`` units.

    Returns
    -----------
    mass_dict : `dict`
        Dictionary containing arrays or floats containing 10-base logarithm
        of halo mass and stellar mass of central galaxies in units of
        either ``h = 1`` or ``h = 0.7``. The units of ``h`` depend on the
        choice of ``return_h0``.

    Notes
    ----------
    The parameter values in Moster+10 were fit to data assuming ``h=0.7``.
    Thus, we will transform our input stellar mass to ``h=0.7`` units,
    evaluate using the Moster parameters, and then transform back to
    ``h=1`` units before returning the result.
    """
    # file_msg = fd.Program_Msg(__file__)
    file_msg = '>>>'
    little_h = 0.72
    ## Checking input parameters
    # `log_halo_mass`
    mstar_valid_types = (int, float, np.ndarray, list)
    if not (isinstance(log_halo_mass, mstar_valid_types)):
        msg = '{0} `log_halo_mass` ({1}) is not a valid type!'.format(
            file_msg, type(log_halo_mass))
        raise LSSUtils_Error(msg)
    ##
    ## Moster dictionary
    param_dict = _retrieve_Moster_default_dict()
    # Converting halo mass from h=1 to h=0.7
    log_halo_mass -= np.log10(little_h)
    halo_mass = (10.**log_halo_mass)
    ##
    ## Moster function
    mass_over_m1 = halo_mass / (10.**param_dict['smhm_m1'])
    m_m1_beta    = mass_over_m1**(-1. * param_dict['smhm_beta'])
    m_m1_gamma   = mass_over_m1**(param_dict['smhm_gamma'])
    inverse_term = (m_m1_beta + m_m1_gamma)**(-1)
    mass_over_M  = 2. * param_dict['smhm_m_M_0'] * inverse_term
    stellar_mass = halo_mass * mass_over_M
    # Converting to log-scale
    log_stellar_mass = np.log10(stellar_mass)
    # Choosing which values to return
    if (return_h0):
        # Returns both masses in units of ``h = 1``
        log_stellar_mass += np.log10(little_h**2)
        log_halo_mass    += np.log10(little_h)
    # Dictionary of masses
    return_obj = [log_halo_mass, log_stellar_mass]
    mass_dict = dict(zip(['m_halo', 'm_stellar'], return_obj))

    return mass_dict

def Moster_example(logmhalo_min=11, logmhalo_max=15, return_h0=True):
    """
    Computes the ``Moster et al. (2010)`` relation and returns a dictionary
    of stellar masses and halo masses of central galaxies

    Parameters
    ------------
    

    return_h0 : `bool`, optional
        If `True`, the function returns both stellar mass and halo mass
        arrays in units of ``h = 1`` units.

    Returns
    ------------
    mass_dict : `dict`
        Dictionary containing arrays or floats containing 10-base logarithm
        of halo mass and stellar mass of central galaxies in units of
        either ``h = 1`` or ``h = 0.7``. The units of ``h`` depend on the
        choice of ``return_h0``.
    """
    # file_msg = fd.Program_Msg(__file__)
    file_msg = '>> '
    # Constants
    dlogM = 1000
    # Checking input parameters
    # `logmhalo_min`
    mhalo_min_max_type = (float, int)
    if not (isinstance(logmhalo_min, mhalo_min_max_type)):
        msg = '{0} `logmhalo_min` ({1}) is not a valid input type: {2}'
        msg = msg.format(file_msg, type(logmhalo_min), mhalo_min_max_type)
        raise TypeError(msg)
    else:
        logmhalo_min = float(logmhalo_min)
    # `logmhalo_max`
    if not (isinstance(logmhalo_max, mhalo_min_max_type)):
        msg = '{0} `logmhalo_max` ({1}) is not a valid input type: {2}'
        msg = msg.format(file_msg, type(logmhalo_max), mhalo_min_max_type)
        raise TypeError(msg)
    else:
        logmhalo_max = float(logmhalo_max)
    # Checking relation between `logmhalo_min` and `mhalo_max`
    if not (logmhalo_min < logmhalo_max):
        msg = '{0} `logmhalo_min` ({1}) must be smaller than `logmhalo_max` '
        msg = '({2})!'
        msg = msg.format(file_msg, logmhalo_min, logmhalo_max)
        raise ValueError(msg)
    # Creating array of halo masses in units of ``h = 1``
    log_halo_mass_arr = np.linspace(logmhalo_min, logmhalo_max, dlogM)
    # Computing stellar masses
    mass_dict = Moster2010_relation(log_halo_mass_arr, return_h0=return_h0)

    return mass_dict
