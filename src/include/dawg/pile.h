#pragma once
#ifndef DAWG_PILE_H
#define DAWG_PILE_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <utility>
#include <string>


#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace dawg{

namespace qi = boost::spirit::qi;
namespace standard = boost::spirit::standard;
namespace phoenix = boost::phoenix;

typedef std::pair<int, int> line_type;

namespace pile {
	typedef std::vector<line_type> subsection_body_type;
	typedef std::pair<std::string, subsection_body_type> subsection_type;
	typedef std::vector<subsection_type> section_body_type;
	typedef std::pair<std::string, std::string> section_header_type;
	typedef std::pair<section_header_type, section_body_type> section_type;
	typedef std::vector<section_type> pile_type;
}


template <typename Iterator>
struct pile_grammar : qi::grammar<Iterator, std::string(), standard::space_type> {
	
	pile_grammar() : pile_grammar::base_type(pile_string) {
		using standard::space; using standard::alnum;
		using standard::graph; using standard::print;
		using boost::spirit::char_;
		using boost::spirit::int_;
		using boost::spirit::raw;
		using boost::spirit::lexeme;
		
		/*
		start = +section;
		section = section_header || section_body;
		section_header = "[[" >> id >> -('=' >> id) >> "]]";
		section_body = +subsection;
		subsection = subsection_header || subsection_body;
		subsection_header = '[' >> id >> ']';
		subsection_body = +line;
		*/
		//line = int_ >> int_; //pile_string; // '"'
		id = +char_("a-zA-Z0-9._") // doesn't work;
		pile_string     %= qqquoted_string |quoted_string | tree_string | bare_string;
		bare_string     %= lexeme[+graph - '"'];
		tree_string     %= lexeme[char_("(") >> +(char_ - ';') >> char_(";")];
		quoted_string   %= lexeme['"' >> *(print - '"') >> '"'];
		qqquoted_string %= lexeme["\"\"\"" >> *(char_ - "\"\"\"") >> "\"\"\""];
		
	}
	
	/*qi::rule<Iterator, std::string(), standard::space_type> start;
	qi::rule<Iterator, pile::section_type(), standard::space_type> section;
	qi::rule<Iterator, pile::section_header_type(), standard::space_type> section_header;
	qi::rule<Iterator, pile::section_body_type(), standard::space_type> section_body;
	qi::rule<Iterator, pile::subsection_type(), standard::space_type> subsection;
	qi::rule<Iterator, std::string(), standard::space_type> subsection_header;
	qi::rule<Iterator, pile::subsection_body_type(), standard::space_type> subsection_body;
	*/
	//qi::rule<Iterator, void(), standard::space_type> line;
	qi::rule<Iterator, std::string(), standard::space_type> id;
	qi::rule<Iterator, std::string(), standard::space_type> pile_string;
	qi::rule<Iterator, std::string(), standard::space_type> bare_string;
	qi::rule<Iterator, std::string(), standard::space_type> tree_string;
	qi::rule<Iterator, std::string(), standard::space_type> quoted_string;
	qi::rule<Iterator, std::string(), standard::space_type> qqquoted_string;	
};


} // namespace dawg

#endif //DAWG_PILE_H

