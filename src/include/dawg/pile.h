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
#include <boost/algorithm/string/predicate.hpp>

#include <utility>
#include <string>
#include <vector>
#include <map>

#include <dawg/log.h>

namespace dawg{

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
	typedef std::vector<section_type> pile_raw_type;
}

struct pile {
	typedef details::pile_raw_type raw;
	struct section {
		typedef std::vector<std::string> value_type;
		typedef std::map<std::string, value_type> map_type;
		std::string inherits;
		map_type db;
	};
	typedef std::map<std::string, section> map_type;
	map_type data;
	
	template<typename Iterator>
	bool parse(Iterator first, Iterator last);
};

template<typename Iterator>
struct white_space : qi::grammar<Iterator> {
	white_space() : white_space::base_type(start) {
		using standard::space; using standard::char_;
		using qi::eol;
		start = space | ('#' >> *(char_ - eol) >> eol);
	}
	qi::rule<Iterator> start;
};

template <typename Iterator, typename skip_type>
struct pile_grammar : qi::grammar<Iterator, pile::raw(), skip_type> {
	pile_grammar() : pile_grammar::base_type(start) {
		using standard::space; using standard::alnum;
		using standard::graph; using standard::print;
		using standard::char_;
		using qi::raw; using qi::lexeme; using qi::eol;
		
		start = +section;
		section = section_header || section_body;
		section_header = "[[" >> -(id >> -('=' >> id)) >> "]]";
		section_body = +subsection;
		subsection = subsection_header || subsection_body;
		subsection_header = '[' >> -id >> ']';
		subsection_body = +line;
		line = id >> '=' >> (pile_string % ',');
		id = lexeme[+(alnum | char_("._-"))];
		pile_string     = qqquoted_string |quoted_string | tree_string | bare_string;
		bare_string     = lexeme[+(graph - char_(",#\"[]=();"))];
		tree_string     = lexeme[char_("(") >> +(char_ - ';') >> char_(";")];
		quoted_string   = lexeme['"' >> *(print - '"') >> '"'];
		qqquoted_string = lexeme["\"\"\"" >> *(char_ - "\"\"\"") >> "\"\"\""];		
	}
	
	qi::rule<Iterator, details::pile_raw_type(), skip_type> start;
	qi::rule<Iterator, details::section_type(), skip_type> section;
	qi::rule<Iterator, details::section_header_type(), skip_type> section_header;
	qi::rule<Iterator, details::section_body_type(), skip_type> section_body;
	qi::rule<Iterator, details::subsection_type(), skip_type> subsection;
	qi::rule<Iterator, details::subsection_body_type(), skip_type> subsection_body;
	qi::rule<Iterator, details::line_type(), skip_type> line;
	qi::rule<Iterator, std::string(), skip_type> subsection_header;
	qi::rule<Iterator, std::string(), skip_type> id;
	qi::rule<Iterator, std::string(), skip_type> pile_string;
	qi::rule<Iterator, std::string(), skip_type> bare_string;
	qi::rule<Iterator, std::string(), skip_type> tree_string;
	qi::rule<Iterator, std::string(), skip_type> quoted_string;
	qi::rule<Iterator, std::string(), skip_type> qqquoted_string;		
};

template<typename Iterator>
bool pile::parse(Iterator first, Iterator last) {
	using boost::algorithm::starts_with;
	raw pyle;
	white_space<Iterator> wasp;
	pile_grammar<Iterator, dawg::white_space<Iterator> > pyler;
	bool r = qi::phrase_parse(first, last, pyler, wasp, pyle);
	if( first != last || !r )
		return DAWG_ERROR("parsing failed.");
	std::string header("_initial_"), subheader("");
	data.clear();
	map_type::iterator sit = data.insert(std::make_pair(header, section())).first;
	sit->second.inherits = "_default_";
	std::pair<map_type::iterator,bool> ret;
	for(raw::const_iterator it = pyle.begin(); it != pyle.end(); ++it) {
		const details::section_type &sec = *it;
		// if section header is not blank, we have a new section
		if(!sec.first.first.empty()) {
			ret = data.insert(std::make_pair(sec.first.first, section()));
			if(!ret.second)
				return DAWG_ERROR("section named " << sec.first.first << " already exists.");
			sit = ret.first;
			// if section parent is black, inherit from previous
			sit->second.inherits = sec.first.second.empty() ? header : sec.first.second;
			// set new header and reset subheader
			header = sec.first.first;
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
				sit->second.db[hh] = lit->second;
			}
		}
			
	}
	
	return true;
}

} // namespace dawg

#endif //DAWG_PILE_H

