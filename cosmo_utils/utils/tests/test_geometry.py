#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018-05-03
# Last Modified: 2018-05-03
from __future__ import print_function, division, absolute_import
__author__     =['Victor Calderon']
__copyright__  =["Copyright 2018 Victor Calderon"]
__email__      =['victor.calderon@vanderbilt.edu']
__maintainer__ =['Victor Calderon']
"""
Set of test functions for the `geometry` functions
"""

## Import modules
import numpy as np
import pytest
from   cosmo_utils.utils import geometry
from   cosmo_utils.utils import file_utils as fd
from   cosmo_utils.custom_exceptions import LSSUtils_Error


## Setting seed
seed = np.random.seed(0)

## Functions

## Testing function `flip_angles`
flip_angles_test_arr = [    (-50, 310.0),
                            (130.5, 130.5),
                            (-0.5, 359.5),
                            ([-10., -20.5, 60.], np.array([350.0, 339.5, 60.0])),
                            (np.array([12., -10]), np.array([12.0, 350.0]))]
@pytest.mark.parametrize('input_ang, output_ang', flip_angles_test_arr)
def test_flip_angles(input_ang, output_ang):
    """
    Tests the function `cosmo_utils.utils.geometry.flip_angles` for input and 
    output parameters

    Parameters
    -----------
    input_ang : float, int, list, `numpy.ndarray`
        Input argument to evaluate

    output_ang : float, int, list, `numpy.ndarray`
        Result to test against. It is the expected result from `flip_angles`
        function.
    """
    ## Checking that outputs from function are the ones to be expected
    output_func = geometry.flip_angles(input_ang)
    # Checking agains 
    np.testing.assert_allclose(output_func, output_ang)

## Testing function Ang_Distance
def test_Ang_Distance_comparison_astropy():
    """
    Tests the function `cosmo_utils.utils.geometry.flip_angles` for input and 
    output parameters.

    This method compares the values obtained from the `haversine` to those 
    using the `astropy` method.
    """
    ## Producing set of Right Ascension and Declination arrays
    ra_lim  = (0, 360.)
    dec_lim = (-90, 90.)
    for ii in range(1, 1000):
        ## Random array for RA and Dec
        ra1  = ra_lim[0]  + np.random.random_sample(ii) * (ra_lim [1]-ra_lim [0])
        dec1 = dec_lim[0] + np.random.random_sample(ii) * (dec_lim[1]-dec_lim[0])
        ra2  = ra_lim[0]  + np.random.random_sample(ii) * (ra_lim [1]-ra_lim [0])
        dec2 = dec_lim[0] + np.random.random_sample(ii) * (dec_lim[1]-dec_lim[0])
        ## Testing outputs from different methods
        out_haversine = geometry.Ang_Distance(ra1, ra2, dec1, dec2,
            method='haversine', unit='deg')
        out_astropy   = geometry.Ang_Distance(ra1, ra2, dec1, dec2,
            method='astropy', unit='deg')
        ## Checking that arrays are the same or similar
        np.testing.assert_allclose(out_haversine, out_astropy)

## Testing `Ang_Distance` for errors - Units
Ang_Distance_test_unit_arr   = [ 'deg2' , 'nan'      ]
@pytest.mark.parametrize('unit', Ang_Distance_test_unit_arr)
def test_Ang_Distance_unit_errors(unit):
    """
    Tests the function `cosmo_utils.utils.geometry.flip_angles` for input and 
    output parameters.

    This function makes sure that errors are raised whenever a wrong 
    input is given.

    Parameters
    -----------
    unit : {'dec','rad'} str, optional
        Unit of `ra1`, `ra2`, `dec1`, and `dec2`.
        This will also determine the final unit that outputs this function.

    method : {'haversine', 'astropy'} str, optional
        Method to use in order to calculate angular separation.
        This variable is to by default to the `haversine` method.
        If `astropy`, it will use the astropy framework to determine the 
        angular separation.
    """
    ## Producing set of Right Ascension and Declination arrays
    ra_lim  = (  0, 360.)
    dec_lim = (-90,  90.)
    for ii in range(1, 1000):
        ## Random array for RA and Dec
        ra1  = ra_lim [0] + np.random.random_sample(ii) * (ra_lim [1]-ra_lim [0])
        dec1 = dec_lim[0] + np.random.random_sample(ii) * (dec_lim[1]-dec_lim[0])
        ra2  = ra_lim [0] + np.random.random_sample(ii) * (ra_lim [1]-ra_lim [0])
        dec2 = dec_lim[0] + np.random.random_sample(ii) * (dec_lim[1]-dec_lim[0])
        ## Testing outputs from different methods
        # Haversine method
        with pytest.raises(LSSUtils_Error):
            out_astropy   = geometry.Ang_Distance(ra1, ra2, dec1, dec2,
                method='astropy', unit=unit)

## Testing `Ang_Distance` for errors - Method
Ang_Distance_test_method_arr = [ 'meter', 'NotMethod']
@pytest.mark.parametrize('method', Ang_Distance_test_method_arr)
def test_Ang_Distance_method_errors(method):
    """
    Tests the function `cosmo_utils.utils.geometry.flip_angles` for input and 
    output parameters.

    This function makes sure that errors are raised whenever a wrong 
    input is given.

    Parameters
    -----------
    unit : {'dec','rad'} str, optional
        Unit of `ra1`, `ra2`, `dec1`, and `dec2`.
        This will also determine the final unit that outputs this function.

    method : {'haversine', 'astropy'} str, optional
        Method to use in order to calculate angular separation.
        This variable is to by default to the `haversine` method.
        If `astropy`, it will use the astropy framework to determine the 
        angular separation.
    """
    ## Producing set of Right Ascension and Declination arrays
    ra_lim  = (  0, 360.)
    dec_lim = (-90,  90.)
    for ii in range(1, 1000):
        ## Random array for RA and Dec
        ra1  = ra_lim [0] + np.random.random_sample(ii) * (ra_lim [1]-ra_lim [0])
        dec1 = dec_lim[0] + np.random.random_sample(ii) * (dec_lim[1]-dec_lim[0])
        ra2  = ra_lim [0] + np.random.random_sample(ii) * (ra_lim [1]-ra_lim [0])
        dec2 = dec_lim[0] + np.random.random_sample(ii) * (dec_lim[1]-dec_lim[0])
        ## Testing outputs from different methods
        # Haversine method
        with pytest.raises(LSSUtils_Error):
            out_astropy   = geometry.Ang_Distance(ra1, ra2, dec1, dec2,
                method=method, unit='deg')



