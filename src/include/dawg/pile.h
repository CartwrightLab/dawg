#pragma once
#ifndef DAWG_PILE_H
#define DAWG_PILE_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/spirit/include/version.hpp>
#if SPIRIT_VERSION < 0x2010
#	error Spirit version 2.1 or greater required.
#endif

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <utility>
#include <string>


namespace dawg{

namespace qi = boost::spirit::qi;
namespace standard = boost::spirit::standard;
namespace phoenix = boost::phoenix;

namespace pile {
	typedef std::pair<std::string, std::vector<std::string> > line_type;
	typedef std::vector<line_type> subsection_body_type;
	typedef std::pair<std::string, subsection_body_type> subsection_type;
	typedef std::vector<subsection_type> section_body_type;
	typedef std::pair<std::string, std::string> section_header_type;
	typedef std::pair<section_header_type, section_body_type> section_type;
	typedef std::vector<section_type> pile_type;
}


template <typename Iterator>
struct pile_grammar : qi::grammar<Iterator, pile::pile_type(), standard::space_type> {
	
	pile_grammar() : pile_grammar::base_type(start) {
		using standard::space; using standard::alnum;
		using standard::graph; using standard::print;
		using standard::char_;
		using qi::raw;
		using qi::lexeme;
		
		start = +section;
		section = section_header || section_body;
		section_header = "[[" >> -(id >> -('=' >> id)) >> "]]";
		section_body = +subsection;
		subsection = subsection_header || subsection_body;
		subsection_header = '[' >> -id >> ']';
		subsection_body = +line;
		line = id >> '=' >> (pile_string % ',');
		id = +(alnum | char_("._-"));
		pile_string     = qqquoted_string |quoted_string | tree_string | bare_string;
		bare_string     = lexeme[+graph - '"'];
		tree_string     = lexeme[char_("(") >> +(char_ - ';') >> char_(";")];
		quoted_string   = lexeme['"' >> *(print - '"') >> '"'];
		qqquoted_string = lexeme["\"\"\"" >> *(char_ - "\"\"\"") >> "\"\"\""];
		
	}
	
	qi::rule<Iterator, pile::pile_type(), standard::space_type> start;
	qi::rule<Iterator, pile::section_type(), standard::space_type> section;
	qi::rule<Iterator, pile::section_header_type(), standard::space_type> section_header;
	qi::rule<Iterator, pile::section_body_type(), standard::space_type> section_body;
	qi::rule<Iterator, pile::subsection_type(), standard::space_type> subsection;
	qi::rule<Iterator, std::string(), standard::space_type> subsection_header;
	qi::rule<Iterator, pile::subsection_body_type(), standard::space_type> subsection_body;
	qi::rule<Iterator, pile::line_type(), standard::space_type> line;
	qi::rule<Iterator, std::string(), standard::space_type> id;
	qi::rule<Iterator, std::string(), standard::space_type> pile_string;
	qi::rule<Iterator, std::string(), standard::space_type> bare_string;
	qi::rule<Iterator, std::string(), standard::space_type> tree_string;
	qi::rule<Iterator, std::string(), standard::space_type> quoted_string;
	qi::rule<Iterator, std::string(), standard::space_type> qqquoted_string;	
};


} // namespace dawg

#endif //DAWG_PILE_H

