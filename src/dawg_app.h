/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 *                                                                          *
 *  This program is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 ****************************************************************************/

#ifndef DAWG_APP_H
#define DAWG_APP_H

#include <vector>
#include <utility>

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace boost {
void validate(boost::any& v, const std::vector<std::string>& xs, boost::tribool*, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
	std::string s(validators::get_single_string(xs, true));

    for (size_t i = 0; i < s.size(); ++i)
        s[i] = char(tolower(s[i]));

    if (s.empty() || s == "on" || s == "yes" || s == "1" || s == "true")
		v = boost::any(boost::tribool(true));
    else if (s == "off" || s == "no" || s == "0" || s == "false")
		v = boost::any(boost::tribool(false));
    else if (s == "null" || s == "maybe" || s == "2" || s == "indeterminate")
		v = boost::any(boost::tribool(boost::indeterminate));
    else
        boost::throw_exception(validation_error(validation_error::invalid_option_value, s));
}
#if !defined(BOOST_NO_STD_WSTRING)
void validate(boost::any& v, const std::vector<std::wstring>& xs, boost::tribool*, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
	std::wstring s(validators::get_single_string(xs, true));

    for (size_t i = 0; i < s.size(); ++i)
        s[i] = char(tolower(s[i]));

    if (s.empty() || s == L"on" || s == L"yes" || s == L"1" || s == L"true")
		v = boost::any(boost::tribool(true));
    else if (s == L"off" || s == L"no" || s == L"0" || s == L"false")
		v = boost::any(boost::tribool(false));
    else if (s == L"null" || s == L"maybe" || s == L"2" || s == L"indeterminate")
		v = boost::any(boost::tribool(boost::indeterminate));
    else
        boost::throw_exception(validation_error(validation_error::invalid_option_value));
}
#endif
}

namespace boost { namespace program_options {
template<>
typed_value<bool>* value(bool* v) {
	return bool_switch(v);
	//typed_value<bool>* r = new typed_value<bool>(v);
    //r->default_value(0, "off");
    //r->implicit_value(1, "on");
	//return r;
}
template<>
typed_value<boost::tribool>* value(boost::tribool* v) {
	//return bool_switch(v);
	typed_value<boost::tribool>* r = new typed_value<boost::tribool>(v);
    r->implicit_value(true, "on");
	return r;
}
}}

/****************************************************************************
 *    class dawg_app                                                        *
 ****************************************************************************/

class dawg_app  {
public:
	dawg_app(int argc, char *argv[]);
	virtual ~dawg_app() { }
	
	virtual int run();

	std::string runname;
	po::options_description desc, indesc;
	po::positional_options_description pdesc;
	po::variables_map vm;

	struct args
	{			
	// use X-Macros to specify argument variables
#	define XM(lname, sname, desc, type, def) type XV(lname) ;
#	include "dawgarg.xmh"
#	undef XM
		std::vector< std::string > input;
	};
	
protected:
	args arg;
};

#endif
