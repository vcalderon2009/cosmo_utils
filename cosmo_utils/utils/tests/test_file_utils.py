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
import os
import numpy as np
import pytest
import cosmo_utils
from cosmo_utils.utils import file_utils as fd
from cosmo_utils.utils import work_paths as wp

## Functions
@pytest.fixture
def get_test_path():
    """
    Path to the Testing directory of `cosmo_utils.utils`

    Returns
    -------
    test_path : str
        Path to the `test` directory
    """
    # Base path
    module_path    = cosmo_utils__path__[0]
    test_path      = os.path.join(module_path,'utils','tests')
    test_data_path = os.path.join(test_path, 'data')
    # test_path      = os.path.join(  base_path,
    #                                 'cosmo_utils',
    #                                 'utils',
    #                                 'tests')
    # test_data_path = os.path.join(test_path, 'data')

    return test_path, test_data_path

##
## Testing the function of Index
test_Index_arr = [  ('.txt', ['data_test.txt']),
                    ('.csv', ['data_test.csv', 'data_test2.csv'])]
@pytest.mark.parametrize('datatype, expected', test_Index_arr)
@pytest.mark.parametr
def test_Index(datatype, expected):
    """
    Tests the function `cosmo_utils.utils.file_utils.Index`

    Parameters
    -----------
    datatype : str
        Type of documents to look for in the directory

    expected : list
        List of files to expect for the given `datatype` in the test path.
    """
    ## Testing different datatypes
    # Main Testing paths
    (   test_path     ,
        test_data_path) = get_test_path()
    # Running command
    out_file  = fd.Index(test_data_path, datatype, basename=True)
    # Checking outputs
    assert(len(out_file) == len(expected))
    assert(np.array_equal(out_file, np.array(expected)))

##
## Testing the function `get_immediate_subdirectories`
def test_get_immediate_subdirectories():
    """
    Tests the function:
        `cosmo_utils.utils.file_utils.get_immediate_subdirectories`
    """
    # Main Testing paths
    (   test_path     ,
        test_data_path) = get_test_path()
    # Running command
    subdir_arr = fd.get_immediate_subdirectories(test_path)
    assert('data' in subdir_arr)

##
## Testing the function `File_Exists`
test_file_exist_arr = [ 'data_test.csv', 'data_test.txt']
@pytest.mark.parametrize('filename', test_file_exist_arr)
def test_File_Exists(filename):
    """
    Tests the function:
        `cosmo_utils.utils.file_utils.File_Exists`

    Parameters
    ----------
    filename : str
        Name of the file being tested for existence
    """
    # Main Testing paths
    (   test_path     ,
        test_data_path) = get_test_path()
    # Running command
    assert(fd.File_Exists(os.path.join(test_data_path, filename)))





