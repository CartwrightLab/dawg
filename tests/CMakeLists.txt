###############################################################################
# Determine the location of the executables/scripts and any data sets required for full package testing
set(DAWG_EXE ${CMAKE_BINARY_DIR}/src/dawg)

###############################################################################
# AddCMakeTest
# A macro to create a test based on a .cmake.in file
# Example: AddCMakeTest(DngCall.Bam "" "Build.dng-call;Build.testdata")
#   - Creates a test called DngCall.Bam
#   - Makes that test depend on Build.dng-call and Build.testdata
#   - Processes the file ./DngCall/BamTest.cmake.in
#   - When the test is called, it calls the .cmake script with possible PREARGS and post ARGs

macro(AddCMakeTest NAME PREARGS DEPENDS)
    string(REPLACE "." "/" pathname "${NAME}")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${pathname}Test.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/${pathname}Test.cmake" @ONLY)
    add_test(NAME ${NAME}
            COMMAND ${CMAKE_COMMAND} ${PREARGS}
            -P "${CMAKE_CURRENT_BINARY_DIR}/${pathname}Test.cmake" ${ARGN})
    if("${DEPENDS}")
        set_tests_properties(${NAME} PROPERTIES DEPENDS "${DEPENDS}")
    endif()
endmacro()

###############################################################################
# Add Tests Here

add_test(Source.CheckBuildDir ${CMAKE_COMMAND} -E md5sum "${CMAKE_SOURCE_DIR}/build/.gitkeep")
set_tests_properties(Source.CheckBuildDir PROPERTIES
        PASS_REGULAR_EXPRESSION "d41d8cd98f00b204e9800998ecf8427e"
        )

add_subdirectory(Unit)
