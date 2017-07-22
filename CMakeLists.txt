###########################################################
#  Dawg's Project Configuration
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0 FATAL_ERROR)
CMAKE_POLICY(SET CMP0028 NEW)
PROJECT(Dawg)

set(NEW_PACKAGE_NAME "dawg")
set(NEW_PACKAGE_VERSION "2-current-rUnknown")

SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Turn on Verbose Makefiles" FORCE)
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)
SET(BUILD_EXTERNAL_PROJECTS OFF CACHE BOOL "Build external dependencies if they are not found.")

SET(USE_STATIC_LIBS OFF CACHE BOOL "Link against static libraries.")
SET(LIBDAWG_USE_STATIC_LIBS OFF CACHE BOOL "Link against a static version of libdawg2.")
SET(RANDOM_GEN "xorshift_64" CACHE STRING "Pseudo-random number generator to use.")
MARK_AS_ADVANCED(RANDOM_GEN)

SET(DAWG_BUILD_DOC OFF CACHE BOOL "TRUE to generate Doxygen API docs, FALSE to ignore it")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

## Guard against bad build-type strings
IF(CMAKE_BUILD_TYPE)
    string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
    if(    NOT cmake_build_type_tolower STREQUAL "debug"
            AND NOT cmake_build_type_tolower STREQUAL "release"
            AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo"
            AND NOT cmake_build_type_tolower STREQUAL "minsizerel"
            AND NOT cmake_build_type_tolower STREQUAL "none")
        message(FATAL_ERROR "Unsupported build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo, MinSizeRel, and None (case-insensitive).")
    endif()
    message(STATUS "Compiling with build type '${CMAKE_BUILD_TYPE}'.")
ENDIF()

#### RPATH Support #####
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

INCLUDE(DawgDevel)
INCLUDE(NewPackage)

if(WIN32 AND NOT UNIX)
	set(DIR_SUFFIX "")
else()
	set(DIR_SUFFIX "\${NEW_PACKAGE_INSTALL_SUFFIX}")
endif()
SET(CMAKE_DATA_DIR "share${DIR_SUFFIX}" CACHE STRING
  "Install location for data (relative to prefix).")
SET(CMAKE_DOC_DIR "doc${DIR_SUFFIX}" CACHE STRING
  "Install location for documentation (relative to prefix).")
SET(CMAKE_MAN_DIR "man" CACHE STRING
  "Install location for man pages (relative to prefix).")
MARK_AS_ADVANCED(CMAKE_DATA_DIR CMAKE_DOC_DIR CMAKE_MAN_DIR)

## Turn on unit and executable testing
ENABLE_TESTING()

################################################################################
# Detect Version Information
FIND_PACKAGE(Git)
SET(GIT_COMMIT_ID "unknown")
IF(GIT_FOUND)
    EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            RESULT_VARIABLE res_var OUTPUT_VARIABLE GIT_COM_ID
            ERROR_QUIET)
    IF(${res_var} EQUAL 0 )
        STRING( REPLACE "\n" "" GIT_COMMIT_ID "${GIT_COM_ID}" )
    ENDIF()
ENDIF()

# try directory name
IF(NOT GIT_COMMIT_ID AND PROJECT_SOURCE_DIR MATCHES "${NEW_PACKAGE_NAME}-([-a-zA-Z0-9.]*)$" )
    SET(GIT_COMMIT_ID "${CMAKE_MATCH_1}")
ENDIF()

# extract project info
IF(GIT_COMMIT_ID MATCHES
        "^v?([0-9]+)\\.([0-9]+)\\.?([0-9]+)?((-[a-zA-Z]+)?(-[0-9]+)?(-g[0-9A-Fa-f]+)?(-dirty)?)$" )
    IF(NOT DEFINED DAWG_VERSION_MAJOR)
        SET(DAWG_VERSION_MAJOR "${CMAKE_MATCH_1}")
    ENDIF()
    IF(NOT DEFINED DAWG_VERSION_MINOR)
        SET(DAWG_VERSION_MINOR "${CMAKE_MATCH_2}")
    ENDIF()
    IF(NOT DEFINED DAWG_VERSION_PATCH)
        SET(DAWG_VERSION_PATCH "${CMAKE_MATCH_3}")
    ENDIF()
    IF(NOT DEFINED DAWG_VERSION_COMMIT)
        STRING(REGEX REPLACE "^-" "" DAWG_VERSION_COMMIT "${CMAKE_MATCH_4}")
    ENDIF()
ELSE()
    # We were unable to determine project information, use the default values
    INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in")
ENDIF()

SET(DAWG_VERSION "${DAWG_VERSION_MAJOR}.${DAWG_VERSION_MINOR}")
SET(DAWG_VERSION_SHORT "${DAWG_VERSION_MAJOR}.${DAWG_VERSION_MINOR}")
IF(DEFINED DAWG_VERSION_PATCH)
    SET(DAWG_VERSION "${DAWG_VERSION}.${DAWG_VERSION_PATCH}")
    SET(DAWG_VERSION_SHORT "${DAWG_VERSION_SHORT}.${DAWG_VERSION_PATCH}")
ENDIF()
IF(DAWG_VERSION_COMMIT)
    SET(DAWG_VERSION "${DAWG_VERSION}-${DAWG_VERSION_COMMIT}")
    SET(DAWG_VERSION_SHORT "${DAWG_VERSION_SHORT}-${DAWG_VERSION_COMMIT}")
ENDIF()

MESSAGE(STATUS "${DAWG_NAME} version: ${DAWG_VERSION}.")

################################################################################
# Packaging Information

SET(CPACK_SOURCE_IGNORE_FILES
  "dawg\\\\.exe"  "dawg$"
)

SET(CPACK_PACKAGE_NAME "Dawg")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Sequence Evolution Simulator")
SET(CPACK_PACKAGE_VENDOR "Reed A. Cartwright")
SET(CPACK_PACKAGE_EXECUTABLES "dawg" "Dawg")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})

IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  SET(CPACK_NSIS_DISPLAY_NAME "Dawg ${NEW_PACKAGE_VERSION}")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\dawg.exe")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\scit.us\\\\projects\\\\dawg\\\\")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\scit.us\\\\projects\\\\dawg\\\\")
  SET(CPACK_NSIS_CONTACT "cartwright@asu.edu")
ELSE()
  SET(CPACK_STRIP_FILES "bin/dawg")
ENDIF()

INCLUDE(CPackBefore)
INCLUDE(CPack)
INCLUDE(ExternalDep)

INSTALL(FILES copying.txt readme.txt changelog.txt DESTINATION ${CMAKE_DOC_DIR})
INSTALL(FILES lambda.pl outsplit.pl nexus2fasta.pl varrep.pl DESTINATION ${CMAKE_DATA_DIR})

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
IF(DAWG_BUILD_DOC) 
    ADD_SUBDIRECTORY(doc)
ENDIF(DAWG_BUILD_DOC)
ADD_SUBDIRECTORY(tests)

ADD_DEPENDENCIES(new_package dawg)
ADD_DEPENDENCIES(dawg version)

########### Add uninstall target ###############
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/Modules/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake") 

#INCLUDE(InstallRequiredSystemLibraries)

