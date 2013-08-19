#pragma once
#ifndef DAWG_TRICK_PARSE_H
#define DAWG_TRICK_PARSE_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <boost/phoenix/core.hpp>

#include <iterator>

#include <dawg/trick.h>
 
#include <boost/spirit/include/version.hpp>
#if SPIRIT_VERSION < 0x2020
#	error Spirit version 2.2 or greater required.
#endif

#include <boost/fusion/include/std_pair.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable: 4127 )
#endif

#include <boost/spirit/include/qi.hpp>

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

namespace dawg {

namespace qi = boost::spirit::qi;
namespace standard = boost::spirit::standard;
namespace phoenix = boost::phoenix;

namespace details {
	typedef std::pair<std::string, std::vector<std::string> > line_type;
	typedef std::vector<line_type> subsection_body_type;
	typedef std::pair<std::string, subsection_body_type> subsection_type;
	typedef std::vector<subsection_type> section_body_type;
	typedef std::pair<std::string, std::string> section_header_type;
	typedef std::pair<section_header_type, section_body_type> section_type;
	typedef std::vector<section_type> trick_raw_type;
}

template<typename Iterator>
struct white_space : qi::grammar<Iterator> {
	white_space() : white_space::base_type(start) {
		using standard::space; using standard::char_;
		using qi::eol; using qi::eoi;
		start = space | ('#' >> *(char_ - eol));
	}
	qi::rule<Iterator> start;
};

template <typename Iterator, typename skip_type>
struct trick_grammar : qi::grammar<Iterator, details::trick_raw_type(), skip_type> {
	trick_grammar() : trick_grammar::base_type(start) {
		using standard::space; using standard::alnum;
		using standard::graph; using standard::print;
		using standard::char_;
		using qi::raw; using qi::lexeme; using qi::eol;
		
		start = *section;
		section = section_header || section_body;
		section_header = "[[" >> id >> -('=' >> id) >> "]]";
		section_body = +subsection;
		subsection = subsection_header || subsection_body;
		subsection_header = '[' >> -id >> ']';
		subsection_body = +line;
		line = id >> '=' >> (trick_string % ',');
		id = lexeme[+(alnum | char_("._-"))];
		trick_string     = qqquoted_string | quoted_string | tree_string | bare_string;
		bare_string     = lexeme[+(graph - char_(",#\"[]=()"))];
		tree_string     = lexeme[char_("(") >> +(char_ - ';') >> char_(";")];
		quoted_string   = lexeme['"' >> *(print - '"') >> '"'];
		qqquoted_string = lexeme["\"\"\"" >> *(char_ - "\"\"\"") >> "\"\"\""];		
	}
	
	qi::rule<Iterator, details::trick_raw_type(), skip_type> start;
	qi::rule<Iterator, details::section_type(), skip_type> section;
	qi::rule<Iterator, details::section_header_type(), skip_type> section_header;
	qi::rule<Iterator, details::section_body_type(), skip_type> section_body;
	qi::rule<Iterator, details::subsection_type(), skip_type> subsection;
	qi::rule<Iterator, details::subsection_body_type(), skip_type> subsection_body;
	qi::rule<Iterator, details::line_type(), skip_type> line;
	qi::rule<Iterator, std::string(), skip_type> subsection_header;
	qi::rule<Iterator, std::string(), skip_type> id;
	qi::rule<Iterator, std::string(), skip_type> trick_string;
	qi::rule<Iterator, std::string(), skip_type> bare_string;
	qi::rule<Iterator, std::string(), skip_type> tree_string;
	qi::rule<Iterator, std::string(), skip_type> quoted_string;
	qi::rule<Iterator, std::string(), skip_type> qqquoted_string;
};

template<typename Iterator>
bool trick::parse(Iterator first, Iterator last) {
	using boost::algorithm::starts_with;
	using boost::algorithm::to_lower;
	details::trick_raw_type pyle;
	white_space<Iterator> wasp;
	trick_grammar<Iterator, dawg::white_space<Iterator> > pyler;
	bool r = qi::phrase_parse(first, last, pyler, wasp, pyle);
	if( first != last || !r )
		return DAWG_ERROR("parsing failed.");
	std::string header("_initial_"), subheader("");
	int autonum = 1; //TODO: What happens with multiple trick files?
	section *psec = &data.front();
	for(details::trick_raw_type::const_iterator it = pyle.begin(); it != pyle.end(); ++it) {
		const details::section_type &sec = *it;
		// if section header is not blank, we have a new section
		if(!sec.first.first.empty()) {
			data.push_back(section());
			psec = &data.back();
			// if section parent is blank, inherit from previous
			psec->inherits = sec.first.second.empty() ? header : sec.first.second;
			// set new header and reset subheader
			if(sec.first.first != "-") {
				header = sec.first.first;
			} else {
				using boost::spirit::karma::generate;
				using boost::spirit::karma::int_;
				char x[16+12] = "Unnamed Section ";
				char *p = x+16;
				generate(p, int_, autonum++);
				*p = 0;
				header = x;
			}
			psec->name = header;
			subheader.clear();
		}
		const details::section_body_type &body = sec.second;
		for(details::section_body_type::const_iterator ssit = body.begin();
			ssit != body.end(); ++ssit) {
			const std::string &h = ssit->first;
			if(h.empty()) // if h is empty, clear subheader
				subheader.clear();
			else if(starts_with(h, "..")) { // if h is '..', strip off tail
				std::string hh = h;
				do { // repeat as needed
					size_t pos = subheader.rfind('.');
					subheader.erase((pos == std::string::npos) ? 0 : pos);
					hh.erase(0, 2);
				} while(starts_with(hh, ".."));
				if(!subheader.empty() && !hh.empty())
					subheader.append(1, '.');
				subheader.append(hh);			
			} else if(starts_with(h, ".")) { // if h begins with a single period, append it to existing subheader
				if(!subheader.empty() && h.size() > 1)
					subheader.append(1, '.').append(h.begin()+1, h.end());
				else
					subheader.append(h);
			} else { // set subheader
				subheader = h;
			}
			// cycle through lines
			for(details::subsection_body_type::const_iterator lit = ssit->second.begin();
				lit != ssit->second.end(); ++lit) {
				std::string hh(subheader);
				if(!hh.empty())
					hh.append(1, '.');
				hh.append(lit->first);
				to_lower(hh);
				psec->db[hh] = lit->second;
			}
		}
	}
	return true;
}

template<typename Char, typename Traits>
inline bool trick::parse_stream(std::basic_istream<Char, Traits>& is) {
	is.unsetf(std::ios::skipws);
	boost::spirit::basic_istream_iterator<Char, Traits> first(is), last;
	return parse(first, last);	
}

} // namespace dawg
 
#endif // DAWG_TRICK_PARSE_H
