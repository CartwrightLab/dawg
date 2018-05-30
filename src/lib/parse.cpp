/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#if _MSC_VER
#	pragma warning(disable: 4127)
#endif

#include <dawg/trick_parse.h>
#include <dawg/wood_parse.h>
#include <dawg/error.h>

#include <iostream>
#include <fstream>

using namespace dawg;

bool trick::parse_file(trick& p, const char *cs) {
	bool ret;
	if(cs == nullptr || strcmp(cs, "")==0 || strcmp(cs, "-")==0) {
		ret = p.parse_stream(std::cin);
	} else {
		std::ifstream is(cs);
		if(!is.is_open()) {
			std::error_code ec = dawg_error::open_input_file_fail;
			DAWG_ERROR_INFO_ = std::string(cs);
			throw ec;
		}
		ret = p.parse_stream(is);
	}
	if(!ret) {
		std::error_code ec = dawg_error::parse_input_file_fail;
		DAWG_ERROR_INFO_ = std::string(cs);
		throw ec;
	}
	return true;
}

bool wood::parse_string(wood &w, const std::string &ss) {
	return w.parse(ss);
}
