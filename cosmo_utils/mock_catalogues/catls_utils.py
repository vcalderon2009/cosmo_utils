#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Victor Calderon
# Created      : 2018-05-08
# Last Modified: 2018-05-08
from __future__ import print_function, division, absolute_import
__author__     =['Victor Calderon']
__copyright__  =["Copyright 2018 Victor Calderon"]
__email__      =['victor.calderon@vanderbilt.edu']
__maintainer__ =['Victor Calderon']
__all__        =[   "catl_keys",
                    "catl_keys_prop"]

## Import modules
import os
import numpy as np
import pandas as pd
from   collections       import Counter
from   cosmo_utils.utils import file_utils   as fd
from   cosmo_utils.utils import work_paths   as wp
from   cosmo_utils.utils import file_readers as fr
from   cosmo_utils.custom_exceptions import LSSUtils_Error

## Functions

## Catalogue Keys - Main
def catl_keys(catl_kind, perf_opt=False, return_type='list'):
    """
    Dictionary keys for the different types of catalogues

    Parameters
    ----------
    catl_kind : {'data', 'mocks'} str, optional
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues

    perf_opt : boolean, optional
        Option for using a `perfect` mock catalogue.

    return_type : {'list', 'dict'} str, optional
        Type of output to the be returned. This variable is set to `list`
        by default.

        Options:
            - 'list' : Returns the values as part of a list
            - 'dict' : Returns the values as part of a python dictionary

    Returns
    ----------
    catl_keys : python dictionary or array_like
        Dictionary/array with the proper keys for the catalogue(s).

        Order : 1) `gm_key`, 2) `id_key`, 3) `galtype_key`

    Raises
    ------------
    LSSUtils_Error : Exception from `LSSUtils_Error`
        Program exception if input parameters are accepted.

    Examples
    ----------
    >>> catl_keys('data', perf_opt=False, return_type='list')
    ['M_h', 'groupid', 'galtype']

    >>> catl_keys('mocks', perf_opt=True, return_type='list')
    ['M_h', 'haloid', 'galtype']

    >>> catl_keys('mocks', perf_opt=True, return_type='dict')
    {'galtype_key': 'galtype', 'gm_key': 'M_h', 'id_key': 'haloid'}
    """
    file_msg = fd.Program_Msg(__file__)
    ## Checking input parameters
    # `catl_kind`
    if not (catl_kind in ['data', 'mocks']):
        msg = '{0} `catl_kind` ({1}) is not a valid input parameter!'.format(
            file_msg, catl_kind)
        raise LSSUtils_Error(msg)
    # `return_type`
    if not (return_type in ['list', 'dict']):
        msg = '{0} `return_type` ({1}) is not a valid input parameter'.format(
            file_msg, return_type)
        raise LSSUtils_Error(msg)
    # `perf_opt`
    if not (isinstance(perf_opt, bool)):
        msg = '{0} `perf_opt` ({1}) must be a boolean object!'.format(
            file_msg, type(perf_opt))
        raise LSSUtils_Error(msg)
    ##
    ## Perfect Catalogue
    if catl_kind == 'data':
        perf_opt = False
    ##
    ## Property keys
    if catl_kind == 'data':
        gm_key, id_key, galtype_key = ['M_h', 'groupid', 'galtype']
    elif catl_kind == 'mocks':
        if perf_opt:
            gm_key, id_key, galtype_key = ['M_h', 'haloid', 'galtype']
        else:
            gm_key, id_key, galtype_key = ['M_group', 'groupid', 'g_galtype']
    ##
    ## Saving values
    if return_type == 'dict':
        catl_objs = {   'gm_key'      : gm_key,
                        'id_key'      : id_key,
                        'galtype_key' : galtype_key}
    elif return_type == 'list':
        catl_objs = [   gm_key, id_key, galtype_key]

    return catl_objs

## Catalogue Keys - Properties
def catl_keys_prop(catl_kind, catl_info='members', return_type='list'):
    """
    Dictionary keys for the diffeent galaxy and group properties of 
    catalogues.

    Parameters
    ------------
    catl_kind : {'data', 'mocks'} str, optional
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues

    catl_info : {'members', 'groups'} str, optional
        Option for which kind of catalogues to use.

        Options:
            - 'members' : Member galaxies of group catalogues
            - 'groups' : Catalogues with `group` information.

    return_type : {'list', 'dict'} str, optional
        Type of output to the be returned. This variable is set to `list`
        by default.

        Options:
            - 'list' : Returns the values as part of a list
            - 'dict' : Returns the values as part of a python dictionary

    Return
    ------------
    catl_objs : python dictionary or array_like
        Dictionary/array with the proper keys for the catalogue(s).

        Order : 1) `ssfr_key`, 2) `mstar_key`

    Raises
    ------------
    LSSUtils_Error : Exception from `LSSUtils_Error`
        Program exception if input parameters are accepted.
    
    Examples
    ------------
    >>> catl_keys_prop('data')
    ['logssfr', 'logMstar_JHU']

    >>> catl_keys_prop('mocks', catl_info='groups', return_type='list')
    ['logssfr', 'logMstar']
    """
    file_msg = fd.Program_Msg(__file__)
    ## Checking input parameters
    catl_kind_valid   = ['data'   , 'mocks' ]
    catl_info_valid   = ['members', 'groups']
    return_type_valid = ['list'   , 'dict'  ]
    # `catl_kind`
    if not (catl_kind in catl_kind_valid):
        msg = '{0} `catl_kind` ({1}) is not a valid input!'.format(
            file_msg, catl_kind)
        raise LSSUtils_Error(msg)
    # `catl_info`
    if not (catl_info in catl_info_valid):
        msg = '{0} `catl_info` ({1}) is not a valid input!'.format(
            file_msg, catl_info)
        raise LSSUtils_Error(msg)
    # `return_type`
    if not (return_type in return_type_valid):
        msg = '{0} `return_type` ({1}) is not a valid input!'.format(
            file_msg, return_type)
        raise LSSUtils_Error(msg)
    ##
    ## Property keys
    ##
    ## Data
    if (catl_kind == 'data'):
        ## Members
        if catl_info == 'members':
            # SSFR and Stellar mass
            logssfr_key, logmstar_key = ['logssfr', 'logMstar_JHU']
        ## Groups
        if catl_info == 'groups':
            # SSFR and Stellar mass
            logssfr_key, logmstar_key = ['logssfr_tot', 'logMstar_tot']
    ##
    ## Mocks
    if (catl_kind == 'mocks'):
        ## Members
        if catl_info == 'members':
            # SSFR and Stellar mass
            logssfr_key, logmstar_key = ['logssfr', 'logMstar']
        ## Groups
        if catl_info == 'groups':
            # SSFR and Stellar mass
            logssfr_key, logmstar_key = ['logssfr', 'logMstar']
    ##
    ## Saving values
    if return_type == 'dict':
        catl_objs = {   'logssfr_key' : logssfr_key ,
                        'logmstar_key': logmstar_key}
    elif return_type == 'list':
        catl_objs = [   logssfr_key, logmstar_key]

    return catl_objs













