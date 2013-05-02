# Portable Release Engineering Script

set(PROJECT_VERSION "") # Default Version
set(PROJECT_NAME "dawg")
set(PROJECT_TITLE "Dawg")
set(PROJECT_DISTS "dawg-2*")
set(SVN_URL "svn://scit.us/${PROJECT_NAME}/")

find_program(SVN_BIN svn)

if(NOT SVN_BIN)
	message(FATAL_ERROR "Could not find Subversion binary.")
endif()

set(ENV{LC_ALL} C)

# Script Options
# RELENG_TAG
# RELENG_M32
# RELENG_M64
# RELENG_TOOLCHAIN

if(NOT RELENG_TAG)
	set(RELENG_TAG "current")
endif()

# Identify Temporary Directory
if(WIN32 AND NOT UNIX)
	set(TMPDIR $ENV{TEMP})
	if(NOT TMPDIR)
		set(TMPDIR "c:/Temp")
	endif()
else()
	set(TMPDIR $ENV{TMPDIR})
	if(NOT TMPDIR)
		set(TMPDIR "/tmp")
	endif()
endif()

string(RANDOM TMP)
set(RELENG_DIR "${TMPDIR}/${PROJECT_NAME}-releng-${TMP}/")
set(ARCHIVE_DIR "${RELENG_DIR}/dawg")

message(STATUS "Using ${RELENG_DIR} to build packages ...")
file(MAKE_DIRECTORY "${RELENG_DIR}")

message(STATUS "Exporting ${PROJECT_TITLE} '${RELENG_TAG}' from SVN...")
execute_process(COMMAND ${SVN_BIN} export "${SVN_URL}/${RELENG_TAG}" "${ARCHIVE_DIR}"
	OUTPUT_VARIABLE SVN_OUTPUT
)

if(SVN_OUTPUT MATCHES "\nExported revision (.*)\\.\n")
	set(PROJECT_SVNREV ${CMAKE_MATCH_1})
else()
	message(FATAL_ERROR "Unable to identify revision")
endif()

message(STATUS "Configuring ${PROJECT_TITLE} '${RELENG_TAG}' ...")
set(CMAKE_DEFS
	-DCMAKE_BUILD_TYPE=Release
	-DUSE_STATIC_LIBS=on
	-DNEW_PACKAGE_VERSION_REV=${PROJECT_SVNREV}
)
set(CMAKE_ARGS "")
if(WIN32 AND NOT UNIX)
	set(CMAKE_ARGS -G "NMake Makefiles" ${CMAKE_ARGS})
	set(MAKE_BIN nmake)
elseif(APPLE)
	SET(CMAKE_DEFS ${CMAKE_DEFS} 
		"-DCMAKE_OSX_ARCHITECTURES=x86_64\\;i386"
		-DCMAKE_OSX_DEPLOYMENT_TARGET=10.5
		-DCPACK_SYSTEM_NAME=Darwin64-universal
	)
	set(MAKE_BIN make)
else()
	set(MAKE_BIN make)	
endif()

if(RELENG_M32)
	set(CMAKE_DEFS ${CMAKE_DEFS} -DCMAKE_C_FLAGS=-m32)
endif()
if(RELENG_M64)
	set(CMAKE_DEFS ${CMAKE_DEFS} -DCMAKE_C_FLAGS=-m64)
endif()
if(RELENG_TOOLCHAIN)
	get_filename_component(RELENG_TOOLCHAIN ${RELENG_TOOLCHAIN} REALPATH)
	set(CMAKE_DEFS ${CMAKE_DEFS} -DCMAKE_TOOLCHAIN_FILE=${RELENG_TOOLCHAIN})
endif()

execute_process(COMMAND ${CMAKE_COMMAND} ${CMAKE_ARGS} .. ${CMAKE_DEFS}
	WORKING_DIRECTORY "${ARCHIVE_DIR}/build")

message(STATUS "Building packages ...")

execute_process(COMMAND ${MAKE_BIN} ${PROJECT_NAME} WORKING_DIRECTORY "${ARCHIVE_DIR}/build")	
execute_process(COMMAND ${MAKE_BIN} new_package WORKING_DIRECTORY "${ARCHIVE_DIR}/build")	
execute_process(COMMAND ${MAKE_BIN} new_package_source WORKING_DIRECTORY "${ARCHIVE_DIR}/build")	

message(STATUS "Relocating packages ...")
file(GLOB DISTS "${ARCHIVE_DIR}/build/${PROJECT_DISTS}")
foreach(dist ${DISTS})
	execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${dist} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

if(NOT NO_CLEAN)
	message(STATUS "Cleaning up ...")
	execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory "${RELENG_DIR}")
endif()

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

