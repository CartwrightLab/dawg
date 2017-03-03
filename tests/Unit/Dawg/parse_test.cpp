#define BOOST_TEST_MODULE Dawg::parse_test
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif // BOOST_TEST_DYN_LINK

#include <string>
#include <array>

#include "../boost_test_helper.h"

#include <dawg/trick.h>
#include <dawg/trick_parse.h>

BOOST_AUTO_TEST_SUITE(test_examples)

static std::string dir_prefix = "tests/Unit/Dawg/";

BOOST_AUTO_TEST_CASE(test_basic_dna)
{
    dawg::trick input;
    std::string exampleFile (dir_prefix + "basic-dna.dawg");
    BOOST_CHECK_EQUAL(dawg::trick::parse_file(input, exampleFile.c_str()), true);
}

BOOST_AUTO_TEST_CASE(test_basic_dna_2)
{
    dawg::trick input;
    std::string exampleFile (dir_prefix + "basic-dna-2.dawg");
    BOOST_CHECK_EQUAL(dawg::trick::parse_file(input, exampleFile.c_str()), true);
}

BOOST_AUTO_TEST_CASE(test_basic_dna_zero_rate)
{
    dawg::trick input;
    std::string exampleFile (dir_prefix + "basic-dna-zero-rate.dawg");
    BOOST_CHECK_EQUAL(dawg::trick::parse_file(input, exampleFile.c_str()), true);
}

BOOST_AUTO_TEST_CASE(test_multiple_models)
{
    dawg::trick input;
    std::string exampleFile (dir_prefix + "multiple-models.dawg");
    BOOST_CHECK_EQUAL(dawg::trick::parse_file(input, exampleFile.c_str()), true);
}

BOOST_AUTO_TEST_CASE(test_recombination)
{
    dawg::trick input;
    std::string exampleFile (dir_prefix + "recombination.dawg");
    BOOST_CHECK_EQUAL(dawg::trick::parse_file(input, exampleFile.c_str()), true);
}

BOOST_AUTO_TEST_CASE(test_segments)
{
    dawg::trick input;
    std::string exampleFile (dir_prefix + "segments.dawg");
    BOOST_CHECK_EQUAL(dawg::trick::parse_file(input, exampleFile.c_str()), true);
}

#if defined(ENABLE_YAML)
BOOST_AUTO_TEST_CASE(test_basic_dna_YAML)
{
    dawg::trick inputDawg, inputYaml;
    std::string dawgExample (dir_prefix + "basic-dna.dawg");
    std::string yamlExample (dir_prefix + "basic-dna.yaml");

    std::array<dawg::trick, 2> inputs;
    inputs.at(0) = inputDawg;
    inputs.at(1) = inputYaml;

    BOOST_CHECK_EQUAL(dawg::trick::parse_file(inputs.at(0), dawgExample.c_str()), true);
    BOOST_CHECK_EQUAL(dawg::trick::parse_file(inputs.at(1), yamlExample.c_str()), true);
    BOOST_CHECK_EQUAL(inputs.at(0).data.size(), inputs.at(1).data.size());

    /* This loop relies on the size assertion above */
    for(unsigned int i = 0; i != inputDawg.data.size(); ++i) {
        BOOST_CHECK_EQUAL(inputDawg.data.at(0).name, inputYaml.data.at(0).name);
        BOOST_CHECK_EQUAL(inputDawg.data.at(0).inherits, inputYaml.data.at(0).inherits);
        /* @TODO : add more robust checks */
    }
}
#endif // defined

BOOST_AUTO_TEST_SUITE_END()