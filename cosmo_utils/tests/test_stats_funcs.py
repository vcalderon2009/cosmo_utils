#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018-04-28
# Last Modified: 2018-04-28
from __future__ import print_function, division, absolute_import
__author__     =['Victor Calderon']
__copyright__  =["Copyright 2018 Victor Calderon"]
__email__      =['victor.calderon@vanderbilt.edu']
__maintainer__ =['Victor Calderon']
"""
Set of test functions for the `stats_func` functions
"""

## Import modules
import numpy as np
import pytest
from   cosmo_utils.utils import stats_funcs

## Functions

##
## Testing `myceil` function
myceil_test = [ (12  , 10, 20.0),
                (13  , 1., 13.0),
                (14.5, 2., 16.0),
                (100 , 10, 100.0),
                (-90 , 10, -90),
                (-89.9, 10, -80),
                (12.05, 0.1, 12.1)]
@pytest.mark.parametrize('x, base, output', myceil_test)
def test_myceil(x, base, output):
    """
    Tests the function `lss_utils.utils.stats_funcs.myceil` for input and 
    output parameters

    Parameters
    ----------
    x : float
        Number to be approximated to closest number to `base`

    base : float
        Base used to calculate the closest largest number

    output : float
        Closest float number to `x`, i.e. upper-bound float
    """
    myceil_out = stats_funcs.myceil(x, base)

    assert( pytest.approx(myceil_out, base) == output)

##
## Testing `myceil` function
myfloor_test = [(12   , 10 , 10.0 ),
                (13   , 1. , 13.0 ),
                (14.5 , 2. , 12.0 ),
                (99   , 20 , 80.0 ),
                (-90  , 10 , -90  ),
                (-89.9, 10 , -90.0),
                (12.05, 0.1, 12.0 )]
@pytest.mark.parametrize('x, base, output', myfloor_test)
def test_myfloor(x, base, output):
    """
    Tests the function `lss_utils.utils.stats_funcs.myceil` for input and 
    output parameters

    Parameters
    ----------
    x : float
        Number to be approximated to closest number to `base`

    base : float
        Base used to calculate the closest largest number

    output : float
        Closest float number to `x`, i.e. upper-bound float
    """
    myceil_out = stats_funcs.myfloor(x, base)

    assert( pytest.approx(myceil_out, base) == output)

##
## Testing `Bins_array_create` function
bins_arr_test = [   ([1,2,3,4,5], 2, [0., 2., 4, 6.]),
                    ([1.2, 2.5, 3, 4], 0.8, [0.8, 1.6, 2.4, 3.2, 4.0]),
                    ([2.25, 3.4, 4.2, 4.3, 2.6], 2, [2., 4, 6.]),
                    ([-6, 2, -8.5, 3.], 3, [-9, -6, -3, 0, 3]),
                    ([1,2,3.5, 4], 2, [0., 2, 4.]),
                    ([-100, 200, 0. -12.5], 50., [-100, -50, 0., 50, 100, 150, 200])]
@pytest.mark.parametrize('arr, base, output', bins_arr_test)
def test_Bins_array_create(arr, base, output):
    """
    Tests the function `lss_utils.utils.stats_funcs.Bins_array_create` for 
    input and output parameters

    Parameters
    -----------
    arr : array_like
        Array of of numbers or floats

    base : int or float, optional
        Interval used to create the evenly-spaced array of elements

    output : `numpy.ndarray`
        Array of elements separated in intervals of `base`
    """
    output = np.array(output)
    # Output from `Bins_array_create` function
    bins_arr = stats_funcs.Bins_array_create(arr, base)
    # Checking that array are equal
    assert(np.allclose(output, bins_arr))

##
## Testing `sigma_calcs`
sigma_calcs_test = [    (1, 100 , (1,100)),
                        (2, 1000, (2,100))]
@pytest.mark.parametrize('return_mean, out_type',[(True, tuple),(False, dict)])
@pytest.mark.parametrize('type_sigma',['std'])
@pytest.mark.parametrize('ndim, nelem, output', sigma_calcs_test)
def test_sigma_calcs_shape(ndim, nelem, output, type_sigma,
    return_mean, out_type):
    """
    Tests the function `lss_utils.utils.stats_funcs.sigma_calcs` for 
    input and output parameters

    Parameters
    -----------
    ndim : int
        Dimesions of the array being tested

    nelem : int
        Number of elements in each sub-array

    output : tuple
        Tuple of elements

    type_sigma : {'std', 'perc'} str
        Type of statistics to use for the different sigmas/St. Dev.

    return_mean : boolean
        If true, it returns `mean`, `std` and `sigmas`.
        If False, it only return `sigmas`.

    out_type : {list, dict} type object
        Type object for the given `return_mean` option
    """
    ## Unpacking input parameters
    ndim_out, len_out = output
    ## Generating dataset
    arr_test = np.random.random((ndim, nelem))
    ## Determining output from `sigma_calcs`
    dict_test      = stats_funcs.sigma_calcs(arr_test, type_sigma=type_sigma,
                        return_mean_std=False)
    dict_test_keys = dict_test.keys()
    ## Testing dimensions and lengths
    for ii in dict_test_keys:
        # Dimensions
        dict_ii_ndim = dict_test[ii].ndim
        assert(dict_ii_ndim == 2)
        # Lengths
        dict_ii_len  = len(dict_test[ii])
        assert(dict_ii_len == 2)
    ##
    ## Testing datatypes
    dict_test_all = stats_funcs.sigma_calcs(arr_test, type_sigma=type_sigma,
                        return_mean_std=return_mean)
    assert(type(dict_test_all) == out_type)
    # Testing types for when `return_mean` == 'True'
    if return_mean:
        assert(type(dict_test_all[0]) ==  dict)
        assert(type(dict_test_all[1]) == np.ndarray)
        assert(type(dict_test_all[2]) == np.ndarray)
















