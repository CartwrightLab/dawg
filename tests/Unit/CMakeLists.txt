set(DAWG_DIR "Dawg")

find_package(TestU01 REQUIRED)
set(TESTU01_LIBS ${TESTU01_LIBRARY})

###############################################################################
# AddUnitTest
# A macro that will add a unit test to the testing queue
# Example: AddUnitTest(dng::seq)
#   - creates the `Unit.dng::seq` test
#   - creates the `unittest_dng_seq` target from the file ./dng/seq.cc
#   - builds the target and runs the right test test-command
#   - the target is linked against Boost::UNIT_TEST_FRAMEWORK and libdng

macro(AddUnitTest NAME)
    string(REPLACE "::" "/" head "${NAME}")
    string(REPLACE "::" "_" exe "${NAME}")
    set(target "unittest_${exe}")
    add_executable("${target}" EXCLUDE_FROM_ALL "${head}.cpp")
    # TODO WIN64: add .exe if needed
    set_target_properties("${target}" PROPERTIES OUTPUT_NAME "${head}")
    get_filename_component(dir "${head}" DIRECTORY)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${dir}")

    target_include_directories("${target}" BEFORE PRIVATE "${CMAKE_SOURCE_DIR}/src/include")

    target_link_libraries("${target}"
            libdawg2
            ${Boost_FILESYSTEM_LIBRARY}
            ${Boost_SYSTEM_LIBRARY}
            ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
            ${TESTU01_LIBS}
            )

    add_test("Unit.${NAME}" ${CMAKE_CTEST_COMMAND}
            --build-and-test "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}"
            --build-generator ${CMAKE_GENERATOR}
            --build-nocmake
            --build-noclean
            --build-project "${CMAKE_PROJECT_NAME}"
            --build-target "${target}"
            --test-command "${CMAKE_CURRENT_BINARY_DIR}/${head}" --log_level=unit_scope
            )
endmacro()

file(COPY "${CMAKE_SOURCE_DIR}/examples/basic-dna.dawg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${DAWG_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/examples/basic-dna-2.dawg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${DAWG_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/examples/basic-dna-zero-rate.dawg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${DAWG_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/examples/multiple-models.dawg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${DAWG_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/examples/pseudogene.dawg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${DAWG_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/examples/recombination.dawg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${DAWG_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/examples/segments.dawg" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${DAWG_DIR}")

###############################################################################
# Add Unit tests based on CXX namespaces

AddUnitTest(${DAWG_DIR}::testu01)
AddUnitTest(${DAWG_DIR}::parse_test)