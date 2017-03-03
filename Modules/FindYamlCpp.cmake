# Locate the yaml-cpp library
#
# This module defines the following variables:
#
# YAMLCPP_LIBRARY      the name of the library;
# YAMLCPP_INCLUDE_DIR  where to find yaml-cpp include files.
# YAMLCPP_FOUND        true if both the YAMLCPP_LIBRARY and YAMLCPP_INCLUDE_DIR have been found.
#
# To help locate the library and include file, you can define a
# variable called YAMLCPP_ROOT which points to the root of the YAMLCPP library
# installation.
#

# default search dirs
set( _yamlcpp_HEADER_SEARCH_DIRS
        "/usr/include"
        "/usr/local/include"
        "C:/Program Files (x86)/yaml-cpp/include" )
set( _yamlcpp_LIB_SEARCH_DIRS
        "/usr/lib"
        "/usr/local/lib"
        "C:/Program Files (x86)/yaml-cpp/lib-msvc110" )

# Check environment for root search directory
set( _yamlcpp_ENV_ROOT $ENV{yampcpp_ROOT} )
if( NOT YAMLCPP_ROOT AND _yamlcpp_ENV_ROOT )
    set(YAMLCPP_ROOT ${_yamlcpp_ENV_ROOT} )
endif()

# Put user specified location at beginning of search
if( YAMLCPP_ROOT )
    list( INSERT _yamlcpp_HEADER_SEARCH_DIRS 0 "${YAMLCPP_ROOT}/include" )
    list( INSERT _yamlcpp_LIB_SEARCH_DIRS 0 "${YAMLCPP_ROOT}/yaml-cpp" )
endif()

# Search for the header
FIND_PATH(YAMLCPP_INCLUDE_DIR "yaml.h"
        PATHS ${_yamlcpp_HEADER_SEARCH_DIRS}
        PATH_SUFFIXES yaml-cpp)

# Search for the library
FIND_LIBRARY(YAMLCPP_LIBRARY NAMES yaml-cpp
        PATHS ${_yamlcpp_LIB_SEARCH_DIRS} )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(YAMLCPP DEFAULT_MSG
        YAMLCPP_LIBRARY YAMLCPP_INCLUDE_DIR)