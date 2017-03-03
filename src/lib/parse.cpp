/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#if _MSC_VER
#	pragma warning(disable: 4127)
#endif

#include <dawg/trick.h>
#include <dawg/trick_parse.h>
#include <dawg/wood_parse.h>

#include <iostream>
#include <fstream>
#include <string>

#include <boost/algorithm/string.hpp>

using namespace dawg;

bool trick::parse_file(trick& p, const char *cs) {
	bool ret;
	if(cs == nullptr || strcmp(cs, "")==0 || strcmp(cs, "-")==0) {
		ret = p.parse_stream(std::cin);
	} else {
        const std::string filename (cs);
        auto pos = filename.find_last_of(".");
        std::string suffix (filename.begin() + pos, filename.end());
        boost::algorithm::to_lower(suffix);
#if defined(ENABLE_YAML)
        if (suffix == ".yaml" || suffix == ".yml") {
            ret = p.parse_yaml(cs);
        }
#endif // defined
        else {
            std::ifstream is(cs);
            if(!is.is_open())
                return DAWG_ERROR("unable to open input file '" << cs << "'");
            ret = p.parse_stream(is);
        }
    }
	if(!ret)
		return DAWG_ERROR("unable to parse input '" << cs << "'");
	return true;
}

bool wood::parse_string(wood &w, const std::string &ss) {
	return w.parse(ss);
}
