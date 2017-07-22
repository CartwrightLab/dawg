# reference: https://github.com/SFML/SFML/blob/master/doc/CMakeLists.txt
# Configure Doxygen commands

# MacOSX checks
IF(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
	# Add some path to search Doxygen in more directories
	SET(ADDITIONAL_PATHS
		/Developer/Applications/Doxygen.app/Contents/Resources
		/Developer/Applications/Doxygen.app/Contents/MacOS
		$ENV{HOME}/Applications/Doxygen.app/Contents/Resources
		$ENV{HOME}/Applications/Doxygen.app/Contents/MacOS
		$ENV{HOME}/Applications/Developer/Doxygen.app/Contents/Resources
		$ENV{HOME}/Applications/Developer/Doxygen.app/Contents/MacOS)

	SET(CMAKE_PROGRAM_PATH ${CMAKE_PROGRAM_PATH} ${ADDITIONAL_PATHS})
ENDIF()

FIND_PACKAGE(Doxygen REQUIRED)

SET(DOXYGEN_INPUT_DIR "${PROJECT_SOURCE_DIR}")
SET(DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/doc")
MESSAGE(STATUS "Doxygen input dir: ${DOXYGEN_INPUT_DIR}")
MESSAGE(STATUS "Doxygen output dir: ${DOXYGEN_OUTPUT_DIR}")

# Windows checks
IF(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
	# If HHC is found, we can generated the CHM (compressed HTML) output
	FIND_PROGRAM(DOXYGEN_HHC_PROGRAM
		NAMES hhc.exe
		PATHS "C:/Program Files/HTML Help Workshop" "C:/Program Files (x86)/HTML Help Workshop"
		DOC "HTML Help Compiler Program")
	IF(DOXYGEN_HHC_PROGRAM)
		SET(DOXYGEN_GENERATE_HTMLHELP YES)
	ELSE()
		SET(DOXYGEN_GENERATE_HTMLHELP NO)
	ENDIF(DOXYGEN_HHC_PROGRAM)
ELSE()
	SET(DOXYGEN_HHC_PROGRAM)
	SET(DOXYGEN_GENERATE_HTMLHELP NO)
ENDIF()

SET(DOXYGEN_CONFIGURED_INPUT "${DOXYGEN_OUTPUT_DIR}/doxyfile")
CONFIGURE_FILE(${DOXYGEN_INPUT_DIR}/doc/doxyfile.in ${DOXYGEN_CONFIGURED_INPUT} @ONLY)
CONFIGURE_FILE(${DOXYGEN_INPUT_DIR}/doc/header.html.in ${DOXYGEN_OUTPUT_DIR}/header.html @ONLY)

# Copy the files needed by the documentation
CONFIGURE_FILE(${DOXYGEN_INPUT_DIR}/doc/stylesheet.css ${DOXYGEN_OUTPUT_DIR}/html/stylesheet.css COPYONLY)

# Target setup
ADD_CUSTOM_TARGET(doc ALL
	COMMAND ${CMAKE_COMMAND} -E echo_append "Building API documentation..."
	COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIGURED_INPUT}
	COMMAND ${CMAKE_COMMAND} -E echo "Done."
	WORKING_DIRECTORY ${DOXYGEN_INPUT_DIR})

# Setup install rules
INSTALL(DIRECTORY ${DOXYGEN_OUTPUT_DIR}/html
	DESTINATION ${INSTALL_MISC_DIR}/doc
	COMPONENT doc)

IF(DOXYGEN_HHC_PROGRAM)
	INSTALL(FILES ${DOXYGEN_OUTPUT_DIR}/dawg.chm
			DESTINATION ${INSTALL_MISC_DIR}/doc
			COMPONENT doc)
ENDIF(DOXYGEN_HHC_PROGRAM)

