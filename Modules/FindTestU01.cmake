# Locate the testu01 library
#
# This module defines the following variables:
#
# TESTU01_LIBRARY      the name of the library;
# TESTU01_INCLUDE_DIR  where to find testu01 include files.
# TESTU01_FOUND        true if both the TESTU01_LIBRARY and TESTU01_INCLUDE_DIR have been found.
#
# To help locate the library and include file, you can define a
# variable called TESTU01_ROOT which points to the root of the TESTU01 library
# installation.
#

# default search dirs
set( _testu01_HEADER_SEARCH_DIRS
        "/usr/include"
        "/usr/local/include"
        "C:/Program Files (x86)/testu01/include" )
set( _testu01_LIB_SEARCH_DIRS
        "/usr/lib"
        "/usr/local/lib"
        "C:/Program Files (x86)/testu01/lib-msvc110" )

# Check environment for root search directory
set( _testu01_ENV_ROOT $ENV{testu01_ROOT} )
if( NOT TESTU01_ROOT AND _testu01_ENV_ROOT )
    set(TESTU01_ROOT ${_testu01_ENV_ROOT} )
endif()

# Put user specified location at beginning of search
if( TESTU01_ROOT )
    list( INSERT _testu01_HEADER_SEARCH_DIRS 0 "${TESTU01_ROOT}/include" )
    list( INSERT _testu01_LIB_SEARCH_DIRS 0 "${TESTU01_ROOT}/testu01" )
endif()

# Search for the header
FIND_PATH(TESTU01_INCLUDE_DIR "TestU01.h"
        PATHS ${_testu01_HEADER_SEARCH_DIRS}
        PATH_SUFFIXES TestU01)

# Search for the library
FIND_LIBRARY(TESTU01_LIBRARY NAMES testu01
        PATHS ${_testu01_LIB_SEARCH_DIRS} )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TESTU01 DEFAULT_MSG
        TESTU01_LIBRARY TESTU01_INCLUDE_DIR)