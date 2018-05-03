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
                            (130.5, 130.5)
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
    np.testing.assert_allclose(actual, desired)

## Testing function Ang_Distance
