###########################################################
#  Dawg's Project Configuration
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

set(NEW_PACKAGE_NAME "dawg")
set(NEW_PACKAGE_VERSION "2-current-rUnknown")

SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Turn on Verbose Makefiles" FORCE)
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

SET(USE_STATIC_LIBS OFF CACHE BOOL "Link against static libraries.")
SET(LIBDAWG_USE_STATIC_LIBS OFF CACHE BOOL "Link against a static version of libdawg2.")
SET(RANDOM_GEN "xorshift_64" CACHE STRING "Pseudo-random number generator to use.")
MARK_AS_ADVANCED(RANDOM_GEN)

IF(USE_STATIC_LIBS)
  SET(Boost_USE_STATIC_LIBS ON)
  SET(Boost_USE_STATIC_RUNTIME ON)
  SET(LIBDAWG_USE_STATIC_LIBS ON FORCE)
ENDIF(USE_STATIC_LIBS)

PROJECT(Dawg)

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

INSTALL(FILES copying.txt readme.txt changelog.txt DESTINATION ${CMAKE_DOC_DIR})
INSTALL(FILES lambda.pl outsplit.pl nexus2fasta.pl varrep.pl DESTINATION ${CMAKE_DATA_DIR})

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(tests)

ADD_DEPENDENCIES(new_package dawg)
ADD_DEPENDENCIES(dawg version)

#INCLUDE(InstallRequiredSystemLibraries)

