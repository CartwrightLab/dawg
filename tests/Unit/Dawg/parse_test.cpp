#define BOOST_TEST_MODULE Dawg::parse_test
#define BOOST_TEST_DYN_LINK

#include <string>

#include "../boost_test_helper.h"
//#include <boost/program_options.hpp>

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

BOOST_AUTO_TEST_SUITE_END()