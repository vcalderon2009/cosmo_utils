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
    array([310.0, 110.0, 350.0])
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
            ang_final = float(ang_converted)
    else:
        try:
            ang = np.asarray(ang)
            # Converting to degrees, if applicable
            if unit == 'rad':
                ang_converted = np.degrees(ang)
            elif unit == 'deg':
                ang_converted = ang
            # Checking the range of `ang`
            ang_converted = np.asarray([xx if xx > 0 else xx+360 for xx in 
                ang_converted])
            # Converting back to radians, if applicable
            if unit == 'rad':
                ang_final = np.radians(ang_converted)
            elif unit == 'deg':
                ang_final = ang_converted
            # Converting to float
            ang_final = ang_final.astype(float)
        except:
            msg = '{0} `ang` could not be converted!'.format(file_msg)
            raise LSSUtils_Error(msg)

    return ang_final

## Calculates the Angular distance between 2 points
def Ang_Distance(ra1, ra2, dec1, dec2, unit='deg', method='haversine'):
    """
    Calculates angular separation between two sets of points with given
    right ascensions and declinations.

    Taken from: https://en.wikipedia.org/wiki/Haversine_formula

    Parameters
    -----------
    ra1, ra2 : float
        Right Ascension of the 1st and 2nd points. Units in `degrees` 
        by default.

    dec1, dec2 : float
        Declination of the 1st and 2nd points. Units in `degrees` by default.
    
    unit : {'dec','rad'} str, optional
        Unit of `ra1`, `ra2`, `dec1`, and `dec2`.
        This will also determine the final unit that outputs this function.

    method : {'haversine', 'astropy'} str, optional
        Method to use in order to calculate angular separation.
        This variable is to by default to the `haversine` method.
        If `astropy`, it will use the astropy framework to determine the 
        angular separation.

    Returns
    -----------
    ang_sep : float
        Angular separation between 1st and 2nd point.
        In units of `degrees`.

    Examples
    -----------
    >>> Ang_Distance(12.0, 25.0, 10.0, -5.0, unit='deg')
    19.8168

    Notes
    -----------
    A = 90. - `dec2`
    B = 90. - `dec1`
    D = `ra1` - `ra2`
    c = Angle between two points
    """
    file_msg = fd.Program_Msg(__file__)
    ## Checking input arguments
    # Units
    if not ((unit == 'deg') or (unit == 'rad')):
        msg = '{0} `unit` ({1}) is not a valid argument'.format(
            file_msg, unit)
        raise LSSUtils_Error(msg)
    # Method
    if not ((method == 'haversine') or (method == 'astropy')):
        msg = '{0} `method` ({1}) is not a valid argument'.format(
            file_msg, method)
        raise LSSUtils_Error(msg)
    ##
    ## Flipping angles
    ra1 = flip_angles(ra1, unit=unit, )
    ra2 = flip_angles(ra2, unit=unit)
    ##
    ## Haversine Method
    if method == 'haversine':
        A = np.radians(90. - dec1)
        B = np.radians(90. - dec2)
        D = np.radians(ra1 - ra2 )
        # Distance
        ang_sep = (np.sin((A-B)*.5))**2. + np.sin(A)*np.sin(B)*(np.sin(D*.5))**2.
        ang_sep = np.degrees(2 * np.arcsin(ang_sep**0.5))
    ##
    ## Astropy Method
    if method == 'astropy':
        # Imports
        from astropy import coordinates as coord
        from astropy.coordinates import SkyCoord
        from astropy import units as u
        # Converting to `units`
        if unit == 'deg':
            unit_opt = u.degree
        elif unit == 'rad':
            unit_opt = u.radians
        # Calculations
        P1    = SkyCoord(ra=ra1, dec=dec1, unit=(unit_opt, unit_opt))
        P2    = SkyCoord(ra=ra2, dec=dec2, unit=(unit_opt, unit_opt))
        ang_sep = P1.separation(P2)
        # Converting to final units
        if unit == 'deg':
            ang_sep = ang_sep.degrees
        elif unit == 'rad':
            ang_sep = ang_sep.radians

    return ang_sep








