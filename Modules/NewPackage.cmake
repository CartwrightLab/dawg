## New Packaging

if(NOT DEFINED NEW_PACKAGE_NAME)
	set(NEW_PACKAGE_NAME "software")
	message(WARNING "NEW_PACKAGE_NAME has defaulted to 'software'.")
endif()

if(NOT DEFINED NEW_PACKAGE_VERSION)
	set(NEW_PACKAGE_VERSION "0.1")
	message(WARNING "NEW_PACKAGE_VERSION has defaulted to '0.1'.")
endif()

if(NOT DEFINED NEW_PACKAGE_VERSION_FILE)
	set(NEW_PACKAGE_VERSION_FILE "version.h")
endif()

# for now use svnversion
find_program(Subversion_SVNVERSION_EXECUTABLE svnversion
  DOC "subversion 'version number' program")
mark_as_advanced(Subversion_SVNVERSION_EXECUTABLE)

if(Subversion_SVNVERSION_EXECUTABLE)
	set(SVNVERSION_FOUND 1)
else()
	set(SVNVERSION_FOUND 0)
endif()

if(Subversion_SVNVERSION_EXECUTABLE)
	execute_process(
	     COMMAND "${Subversion_SVNVERSION_EXECUTABLE}" -n "${CMAKE_SOURCE_DIR}"
	     OUTPUT_VARIABLE SVN_OUTPUT
	     ERROR_QUIET
	)
	if(SVN_OUTPUT MATCHES "^[0-9]")
		string(REGEX REPLACE "(-rUnknown|-r[0-9].*)$" "-r${SVN_OUTPUT}"
			NEW_PACKAGE_VERSION "${NEW_PACKAGE_VERSION}")
	endif()
	string(REPLACE ":" "_" NEW_PACKAGE_VERSION "${NEW_PACKAGE_VERSION}")
endif()

################################################################################
# Write helper files to the binary directory

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/newpkg/")

#### version.cmake ####
FILE(WRITE "${CMAKE_BINARY_DIR}/newpkg/version.cmake" "
set(NEW_PACKAGE_VERSION \"${NEW_PACKAGE_VERSION}\")
set(NEW_PACKAGE_NAME \"${NEW_PACKAGE_NAME}\")

if(${SVNVERSION_FOUND})
	execute_process(
	     COMMAND \"${Subversion_SVNVERSION_EXECUTABLE}\" -n \"${CMAKE_SOURCE_DIR}\"
	     OUTPUT_VARIABLE SVN_OUTPUT
	     ERROR_QUIET
	)
	if(SVN_OUTPUT MATCHES \"^[0-9]\")
		string(REGEX REPLACE \"(-rUnknown|-r[0-9].*)\$\" \"-r\${SVN_OUTPUT}\"
			NEW_PACKAGE_VERSION \"\${NEW_PACKAGE_VERSION}\")
#		set(NEW_PACKAGE_VERSION \"\${NEW_PACKAGE_VERSION}-r\${SVN_OUTPUT}\")
	endif()
	string(REPLACE \":\" \"_\" NEW_PACKAGE_VERSION \"\${NEW_PACKAGE_VERSION}\")
endif()

# Prefer .pkg
if(EXISTS ${CMAKE_SOURCE_DIR}/${NEW_PACKAGE_VERSION_FILE}.newpkg)
	set(INFILE ${CMAKE_SOURCE_DIR}/${NEW_PACKAGE_VERSION_FILE}.newpkg)
else()
	set(INFILE ${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE}.in)
endif()

configure_file(\${INFILE} ${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE} @ONLY)
")
####

#### version.h.in ####
file(WRITE "${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE}.in"
"#pragma once
#ifndef NEW_PACKAGE_VERSION_H
#define NEW_PACKAGE_VERSION_H

#if 0
SET(NEW_PACKAGE_NAME    \"\@NEW_PACKAGE_NAME\@\")
SET(NEW_PACKAGE_VERSION \"\@NEW_PACKAGE_VERSION\@\")
#endif

#define NEW_PACKAGE_NAME    \"\@NEW_PACKAGE_NAME\@\"
#define NEW_PACKAGE_VERSION \"\@NEW_PACKAGE_VERSION\@\"

#define NEW_PACKAGE_STRING (NEW_PACKAGE_NAME \" \" NEW_PACKAGE_VERSION)

#endif
")
####

#### new_install.cmake ####
file(WRITE "${CMAKE_BINARY_DIR}/newpkg/install.cmake" "
if(NOT CPACK_INSTALL_CMAKE_PROJECTS)
	if(EXISTS \"${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE}\")
		execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
			\"${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE}\"
	        \"\${CMAKE_INSTALL_PREFIX}/${NEW_PACKAGE_VERSION_FILE}.newpkg\"
		)
	else()
		set(NEW_PACKAGE_NAME \"\${CPACK_PACKAGE_NAME}\")
		set(NEW_PACKAGE_VERSION \"\${CPACK_PACKAGE_VERSION}\")
		configure_file(\"${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE}.in\"
			\"\${CMAKE_INSTALL_PREFIX}/${NEW_PACKAGE_VERSION_FILE}.newpkg\"
			@ONLY)
	endif()
	file(MAKE_DIRECTORY \"\${CMAKE_INSTALL_PREFIX}/build\")
	file(WRITE \"\${CMAKE_INSTALL_PREFIX}/build/.empty\" "")
endif()
")
####

#### new_package.cmake ####
file(WRITE "${CMAKE_BINARY_DIR}/newpkg/new_package.cmake" "
INCLUDE(CPackConfig.cmake)

SET(NEW_PACKAGE_VERSION \"\${CPACK_PACKAGE_VERSION}\")
SET(NEW_PACKAGE_NAME \"\${CPACK_PACKAGE_NAME}\")

INCLUDE(\"${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE}\")

# Some files use directories that contain the version number
# Update as needed
IF(NEW_PACKAGE_VERSION MATCHES \"^([0-9]+)\\\\.\")
	SET(ENV{NEW_PACKAGE_INSTALL_DIR_VERSION} \${CMAKE_MATCH_1})
ENDIF()

EXECUTE_PROCESS(COMMAND \${CMAKE_CPACK_COMMAND}
	-D NEW_PACKAGE=TRUE
	-D CPACK_PACKAGE_FILE_NAME=\${NEW_PACKAGE_NAME}-\${NEW_PACKAGE_VERSION}-\${CPACK_SYSTEM_NAME}
	--config CPackConfig.cmake
	-P \"\${NEW_PACKAGE_NAME}\"
	-R \"\${NEW_PACKAGE_VERSION}\"
	TIMEOUT 3600
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
")
####

#### new_package_source.cmake ####
file(WRITE "${CMAKE_BINARY_DIR}/newpkg/new_package_source.cmake" "
INCLUDE(CPackSourceConfig.cmake)

SET(NEW_PACKAGE_VERSION \"\${CPACK_PACKAGE_VERSION}\")
SET(NEW_PACKAGE_NAME \"\${CPACK_PACKAGE_NAME}\")

INCLUDE(\"${CMAKE_BINARY_DIR}/${NEW_PACKAGE_VERSION_FILE}\")

EXECUTE_PROCESS(COMMAND \${CMAKE_CPACK_COMMAND}
	-D NEW_PACKAGE=TRUE
	-D CPACK_PACKAGE_FILE_NAME=\${NEW_PACKAGE_NAME}-\${NEW_PACKAGE_VERSION}-\${CPACK_SYSTEM_NAME}
	--config CPackSourceConfig.cmake
	-P \"\${NEW_PACKAGE_NAME}\"
	-R \"\${NEW_PACKAGE_VERSION}\"
	TIMEOUT 3600
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
")
####

################################################################################
SET(CPACK_INSTALL_SCRIPT "${CMAKE_BINARY_DIR}/newpkg/install.cmake")

install(CODE "
if(DEFINED ENV{NEW_PACKAGE_INSTALL_DIR_VERSION})
	set(NEW_PACKAGE_INSTALL_DIR_VERSION \"\$ENV{NEW_PACKAGE_INSTALL_DIR_VERSION}\")
else()
	set(NEW_PACKAGE_INSTALL_DIR_VERSION \"${NEW_PACKAGE_INSTALL_DIR_VERSION}\")
endif()
")

if(WIN32 AND NOT UNIX)
  set(NEW_PACKAGE_DIR)
else()
  set(NEW_PACKAGE_DIR "/${NEW_PACKAGE_NAME}-\${NEW_PACKAGE_INSTALL_DIR_VERSION}")
endif()

ADD_CUSTOM_TARGET(version
    COMMAND ${CMAKE_COMMAND}
    	 -P ${CMAKE_BINARY_DIR}/newpkg/version.cmake
)

ADD_CUSTOM_TARGET(new_package
	COMMAND ${CMAKE_COMMAND} 
		 -P ${CMAKE_BINARY_DIR}/newpkg/new_package.cmake
)

ADD_CUSTOM_TARGET(new_package_source
	COMMAND ${CMAKE_COMMAND} 
		 -P ${CMAKE_BINARY_DIR}/newpkg/new_package_source.cmake
)

ADD_DEPENDENCIES(new_package version)
ADD_DEPENDENCIES(new_package_source version)

set(CPACK_PACKAGE_NAME    "${NEW_PACKAGE_NAME}")
set(CPACK_PACKAGE_VERSION "${NEW_PACKAGE_VERSION}")

