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
    {'gm_key': 'M_h', 'id_key': 'haloid', 'galtype_key': 'galtype'}
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
            - `members` : Member galaxies of group catalogues
            - `groups` : Catalogues with `group` information.

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

## Extracting path of synthetic catalogues
def catl_sdss_dir(catl_kind='data', catl_type='mr', sample_s='19',
    catl_info='members', halotype='fof', clf_method=3, hod_n=0, clf_seed=1235,
    perf_opt=False, print_filedir=True):
    """
    Extracts the path to the synthetic catalogues.

    Parameters
    -----------
    catl_kind : {'data', 'mocks'} str, optional
        Type of catalogue to use. This variable is set to `data` by default.

        Options:
            - `data` : catalogues come from SDSS `real` catalogue
            - `mocks` : catalogue come from SDSS `mock` catalogues

    catl_type : {'mr', 'mstar'} str, optional
        Type of catalogue to use. It shows which abundance matching method
        was used for the CLF when assigning halo masses. This variable is 
        set to 'mr' by default.

        Options:
            - `mr` : Uses r-band absolute magnitude
            - `mstar` : Uses stellar masses

    sample_s : {'19', '20', '21'} str, optional
        Volume-limited sample to use. This variable is set to '19' by default.

        Options:
            - '19' : Uses the Mr19 volume-limited sample, i.e. 'Consuelo'
            - '20' : Uses the Mr20 volume-limited sample, i.e. 'Esmeralda'
            - '21' : Uses the Mr21 volume-limited sample, i.e. 'Carmen'

    catl_info : {'members', 'groups'} str, optional
        Option for which kind of catalogues to use.

        Options:
            - `members` : Member galaxies of group catalogues
            - `groups` : Catalogues with `group` information.

    halotype : {'fof', 'so'} str, optional
        Type of the dark matter halo of the simulation used to create the 
        synthetic catalogues. This variable is set to `fof` by default.

        Options:
            - 'fof': Friends-of-Friends halos.
            - 'so' : Spherical overdensity halos.

    clf_method : {1, 2, 3} int, optional
        Method for assigning galaxy properties to mock galaxies.
        This variable is set to `3` by default.

        Options:
            - `1` : Independent assigment of (g-r) color, sersic, and log(ssfr)
            - `2` : (g-r) decides active/passive designation and draw values 
                    independently.
            - `3` : (g-r) decides active/passive designations, and 
                    assigns other galaxy properties for that given galaxy.

    hod_n : {0, 1} int, optional
        HOD model to use. Only relevant when `catl_kind == mocks`.

    clf_seed : int, optional
        Seed used for the `CLF` random seed. This variable is set to `1235` 
        by default.

    perf_opt : boolean, optional
        If True, it chooses to analyze the `perfect` set of synthetic
        catalogues. This variable is set to `False` by default.

    print_filedir : boolean, optional
        If True, the output directory is printed onto the screen.
    
    Returns
    -----------
    catls_path : str
        Path to the desired set of synthetic catalogues.

    Raises
    ------------
    LSSUtils_Error : Exception from `LSSUtils_Error`
        Program exception if input parameters are accepted.
    """
    file_msg = fd.Program_Msg(__file__)
    ## Checking input parameters
    catl_kind_valid  = ['data', 'mocks' ]
    catl_type_valid  = ['mr', 'mstar']
    sample_s_valid   = ['19', '20', '21']
    catl_info_valid  = ['members', 'groups']
    halotype_valid   = ['fof', 'so']
    clf_method_valid = [1, 2, 3]
    hod_n_valid      = [0, 1]
    # `catl_kind`
    if not (catl_kind in catl_kind_valid):
        msg = '{0} `catl_kind` ({1}) is not a valid input!'.format(file_msg,
            catl_kind)
        raise LSSUtils_Error(msg)
    # `catl_type`
    if not (catl_type in catl_type_valid):
        msg = '{0} `catl_type` ({1}) is not a valid input!'.format(file_msg,
            catl_type)
        raise LSSUtils_Error(msg)
    # `sample_s`
    if not (sample_s in sample_s_valid):
        msg = '{0} `sample_s` ({1}) is not a valid input!'.format(file_msg,
            sample_s)
        raise LSSUtils_Error(msg)
    # `catl_info`
    if not (catl_info in catl_info_valid):
        msg = '{0} `catl_info` ({1}) is not a valid input!'.format(file_msg,
            catl_info)
        raise LSSUtils_Error(msg)
    # `halotype`
    if not (halotype in halotype_valid):
        msg = '{0} `halotype` ({1}) is not a valid input!'.format(file_msg,
            halotype)
        raise LSSUtils_Error(msg)
    # `clf_method`
    if not (clf_method in clf_method_valid):
        msg = '{0} `clf_method` ({1}) is not a valid input!'.format(file_msg,
            clf_method)
        raise LSSUtils_Error(msg)
    # `hod_n`
    if not (hod_n in hod_n_valid):
        msg = '{0} `hod_n` ({1}) is not a valid input!'.format(file_msg,
            hod_n)
        raise LSSUtils_Error(msg)
    # `perf_opt`
    if not (isinstance(perf_opt, bool)):
        msg = '{0} `perf_opt` ({1}) is not a valid type!'.format(file_msg,
            type(perf_opt))
        raise LSSUtils_Error(msg)
    # `print_filedir`
    if not (isinstance(print_filedir, bool)):
        msg = '{0} `print_filedir` ({1}) is not a valid type!'.format(file_msg,
            type(print_filedir))
        raise LSSUtils_Error(msg)
    ##
    ## Type of catalogue
    if catl_info == 'members':
        catl_info_str = 'member_galaxy_catalogues'
    elif catl_info == 'groups':
        catl_info_str = 'group_galaxy_catalogues'
    ##
    ## Perfect catalogue
    if perf_opt:
        # Data
        if catl_kind == 'data':
            msg = '{0} Invalid `catl_kind` ({1}) for when `perf_opt == True'
            msg = msg.format(file_msg, catl_kind)
            raise LSSUtils_Error(msg)
        # Mocks
        catl_info_perf_str = 'perfect_{0}'.format(catl_info_str)
    else:
        # Mocks
        catl_info_perf_str = catl_info_str
    ##
    ## Extracting path of the files
    # Data
    if catl_kind == 'data':
        # Joining paths
        filedir = os.path.join( wp.get_output_path(),
                                'SDSS',
                                catl_kind,
                                catl_type,
                                'Mr' + sample_s,
                                catl_info_perf_str)
    # Mocks
    if catl_kind == 'mocks':
        # Joining paths
        filedir = os.path.join( wp.get_output_path(),
                                'SDSS',
                                catl_kind,
                                'halos_{0}'.format(halotype),
                                'hod_model_{0}'.format(hod_n),
                                'clf_seed_{0}'.format(clf_seed),
                                'clf_method_{0}'.format(clf_method),
                                catl_type,
                                'Mr' + sample_s,
                                catl_info_perf_str)
    ##
    ## Making sure `filedir` exists
    if not (os.path.exists(filedir)):
        msg = '{0} `filedir` ({1}) does NOT exist! Check input variables'
        msg = msg.format(file_msg, filedir)
        raise LSSUtils_Error(msg)
    ##
    ## Printing out paths
    if print_filedir:
        print('{0} `filedir`: {1}'.format(file_msg, filedir))

    return filedir
        


















