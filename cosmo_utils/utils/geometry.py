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
__all__        =[   "flip_angles"]
"""
Set of geometrical definitions for translations, coordinate tranformations, 
etc.
"""

## Import modules
import numpy as np
from   cosmo_utils.utils import file_utils as fd
from   cosmo_utils.custom_exceptions import LSSUtils_Error

## Functions

# Restricting angles to be between 0 and 360
def flip_angles(ang, unit='deg'):
    """
    Ensures that an angle is always between 0 and 360 degrees.

    Parameters
    ----------
    ang : float, int, `numpy.ndarray`
        Angle in units of degrees.

    unit : {'deg', 'rad'} str, optional
        Unit of the angle. This variable is set to 'deg' by default.
        If 'rad', the output will be in units of radians.

    Returns
    ----------
    ang_out : float
        Convertions of `ang`, ranging from 0 to 360 degrees.

    Raises
    ----------
    LSSUtils_Error : Exception
        This error gets raised when `ang` is not a digit or a number.

    Examples
    ----------
    >>> flip_angles(-50, unit='deg')
    310.0
    
    >>> flip_angles(110, unit='deg')
    110.0

    >>> flip_angles([-50, 110, -10], unit='deg')
    array([-50.0, 110.0, 350.0])
    """
    file_msg = fd.Program_Msg(__file__)
    # Checking type of `ang`
    if not ((type(ang) == float) or (type(ang) == int) or 
            (type(ang) == list ) or (type(ang) == np.ndarray)):
        msg = '{0} `ang` ({1}) is not a number! Exiting!'.format(
            file_msg, ang)
    # Checking `unit`
    if not (unit.lower() in ['deg', 'rad']):
        msg = '{0} `unit` ({1}) is not a valid option! Exiting!'.format(
            file_msg, unit)
    ##
    ## Converting angle
    if isinstance(ang, float) or isinstance(ang, int):
        ang = float(ang)
        # Converting from radians to degrees, if applicable
        if unit == 'rad':
            ang_converted = np.degrees(ang)
        elif unit == 'deg':
            ang_converted = ang
        # Checking the range of `ang`
        if ang_converted < 0:
            ang_converted += 360
        # Converting back to radians, if applicable
        if unit == 'rad':
            ang_final = np.radians(ang_converted)
        elif unit == 'deg':
            ang_final = ang_converted
    else:
        try:
            ang = np.asarray(ang)
            # Converting to degrees, if applicable
            if unit == 'rad':
                ang_converted = np.degrees(ang)
            elif unit == 'deg':
                ang_converted = ang
            # Checking the range of `ang`
            ang_converted = np.asarray([xx if xx > 0 else x+360 for xx in 
                ang_converted])
            # Converting back to radians, if applicable
            if unit == 'rad':
                ang_final = np.radians(ang_converted)
            elif unit == 'deg':
                ang_final = ang_converted
        except:
            msg = '{0} `ang` could not be converted!'.format(file_msg)
            raise LSSUtils_Error(msg)

    return ang_final








