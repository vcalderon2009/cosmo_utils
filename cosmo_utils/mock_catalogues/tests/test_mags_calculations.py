#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018-05-07
# Last Modified: 2018-05-07
from __future__ import print_function, division, absolute_import
__author__     =['Victor Calderon']
__copyright__  =["Copyright 2018 Victor Calderon"]
__email__      =['victor.calderon@vanderbilt.edu']
__maintainer__ =['Victor Calderon']
"""
Set of test functions for the `mags_calculations` functions
"""

## Import modules
import numpy as np
import pytest
from   cosmo_utils.mock_catalogues import mags_calculations
from   cosmo_utils.custom_exceptions import LSSUtils_Error

## Functions

## Testing `apparent_to_absolute_magnitude` function
app_abs_test_arr = [    (-0.72, -3.1 , 30.1      , 'pc'),
                        ( 0.14, -7.1 , 0.2761    , 'kpc'),
                        ( 1.26, -7.1 , 0.4908    , 'kpc')]
@pytest.mark.parametrize('app_mag, abs_mag, lum_dist, unit', app_abs_test_arr)
def test_apparent_to_absolute_magnitude(app_mag, abs_mag, lum_dist, unit):
    """
    Tests the function 
    `cosmo_utils.mock_catalogues.mags_calculations.apparent_to_absolute_magnitude`
    for input and output parameters.

    This function tests whether or not the function recovers the 
    expected absolute magnitude.

    Parameters
    ----------
    app_mag : array_like
        Array of apparent magnitude(s)

    abs_mag : np.ndarray
        Array of absolute magnitudes. `abs_mag` is a float if 
        `app_mag` is a float or int.

    lum_dist : array_like
        Array of luminosity distnace to object. In units of `Mpc`.

    unit : {'pc', 'kpc', 'mpc'} str, optional
        Unit to use for `lum_dist`. This variable is set to `mpc` by
        default. When `pc`, the units are in parsecs, while `mpc` is for 
        distances in mega-parsecs, etc.
    """
    ## Calculations
    out_abs_mag = mags_calculations.apparent_to_absolute_magnitude(
        app_mag, lum_dist, unit=unit)
    ## Comparing to expected absolute distance
    np.testing.assert_almost_equal(abs_mag, out_abs_mag, decimal=2)

