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
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iterator>

#include <dawg/log.h>
#include <dawg/utils/foreach.h>

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
	typedef std::vector<section_type> pile_raw_type;
}

struct pile {
	typedef details::pile_raw_type raw;
	struct section {
		typedef std::vector<std::string> value_type;
		typedef std::map<std::string, value_type> db_type;
		std::string name;
		std::string inherits;
		db_type db;
		
		template<typename T>
		inline void get(const std::string& k, T& r) const;
		template<typename T, typename A>
		inline void get(const std::string& k, std::vector<T,A>& r) const;
		
		inline void read_aliases();
		
	private:
		static inline void conv(const std::string& ss, std::string& r);
		static inline void conv(const std::string& ss, double& r);
		static inline void conv(const std::string& ss, bool& r);
		static inline void conv(const std::string& ss, unsigned int& r);
		
		inline void read_alias(const std::string& a, const std::string& b);
	};
	typedef std::vector<section> data_type;
	data_type data;

	inline bool parse_file(const char *cs);
	template<typename Iterator>
	bool parse(Iterator first, Iterator last);
	template<typename Char, typename Traits>
	inline bool parse_stream(std::basic_istream<Char, Traits>& is);
	
	pile() {
		data.push_back(section());
		data.back().name = "_initial_";
		data.back().inherits = "_default_";
	}
	
	inline void read_aliases();	
};

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
struct pile_grammar : qi::grammar<Iterator, pile::raw(), skip_type> {
	pile_grammar() : pile_grammar::base_type(start) {
		using standard::space; using standard::alnum;
		using standard::graph; using standard::print;
		using standard::char_;
		using qi::raw; using qi::lexeme; using qi::eol;
		
		start = +section;
		section = section_header || section_body;
		section_header = "[[" >> id >> -('=' >> id) >> "]]";
		section_body = +subsection;
		subsection = subsection_header || subsection_body;
		subsection_header = '[' >> -id >> ']';
		subsection_body = +line;
		line = id >> '=' >> (pile_string % ',');
		id = lexeme[+(alnum | char_("._-"))];
		pile_string     = qqquoted_string | quoted_string | tree_string | bare_string;
		bare_string     = lexeme[+(graph - char_(",#\"[]=()"))];
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
	using boost::algorithm::to_lower;
	raw pyle;
	white_space<Iterator> wasp;
	pile_grammar<Iterator, dawg::white_space<Iterator> > pyler;
	bool r = qi::phrase_parse(first, last, pyler, wasp, pyle);
	if( first != last || !r )
		return DAWG_ERROR("parsing failed.");
	std::string header("_initial_"), subheader("");
	int autonum = 1;
	section *psec = &data.front();
	for(raw::const_iterator it = pyle.begin(); it != pyle.end(); ++it) {
		const details::section_type &sec = *it;
		// if section header is not blank, we have a new section
		if(!sec.first.first.empty()) {
			data.push_back(section());
			psec = &data.back();
			// if section parent is black, inherit from previous
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

// TODO: Convert to newer 2.2 framework?
template<typename Char, typename Traits>
inline bool pile::parse_stream(std::basic_istream<Char, Traits>& is) {
	is.unsetf(std::ios::skipws);
	typedef std::istreambuf_iterator<Char, Traits> base_iterator_type;
	boost::spirit::multi_pass<base_iterator_type> first = 
		boost::spirit::make_default_multi_pass(base_iterator_type(is));
	boost::spirit::multi_pass<base_iterator_type> last = 
		boost::spirit::make_default_multi_pass(base_iterator_type());
	return parse(first, last);	
}

inline bool pile::parse_file(const char *cs) {
	bool ret;
	if(cs == NULL || strcmp(cs, "")==0 || strcmp(cs, "-")==0) {
		ret = parse_stream(std::cin);
	} else {
		std::ifstream is(cs);
		if(!is.is_open())
			return DAWG_ERROR("unable to open input file '" << cs << "'");
		ret = parse_stream(is);
	}
	if(!ret)
		return DAWG_ERROR("unable to parse input '" << cs << "'");
	return true;
}

template<typename T>
inline void pile::section::get(const std::string& k, T& r) const {
	db_type::const_iterator it;
	if((it = db.find(k)) != db.end() && !it->second.empty()) {
		section::conv(it->second.front(), r);
	}
}

template<typename T, typename A>
inline void pile::section::get(const std::string& k, std::vector<T,A>& r) const {
	db_type::const_iterator it;
	if((it = db.find(k)) != db.end()) {
		T x;
		r.clear();
		foreach(const std::string &ss, it->second) {
			section::conv(ss, x);
			r.push_back(x);
		}
	}
}

inline void pile::section::conv(const std::string& ss, std::string& r) {
	r = ss;
}

inline void pile::section::conv(const std::string& ss, double& r) {
	r = strtod(ss.c_str(), NULL);
}

inline void pile::section::conv(const std::string& ss, unsigned int& r) {
	r = strtoul(ss.c_str(), NULL, 0);
}

// A value is false if it is equal to 0, f, false, off, no, or blank
inline void pile::section::conv(const std::string& ss, bool& r) {
	using boost::algorithm::iequals;
	r = !(ss.empty() || iequals(ss, "false") || iequals(ss, "0") || iequals(ss, "f")
		|| iequals(ss, "off") || iequals(ss, "no"));
}

inline void pile::section::read_alias(const std::string& a, const std::string& b) {
	db_type::const_iterator it;
	// if b exists or a doesn't, stop
	if(db.find(b) != db.end() || (it = db.find(a)) == db.end())
		return;
	db[b] = it->second;
}

inline void pile::section::read_aliases() {
	#define XM(aname, bname) read_alias(_P(aname), _P(bname));
	#include <dawg/details/aliases.xmh>
	#undef XM
}

inline void pile::read_aliases() {
	foreach(section &sec, data) {
		sec.read_aliases();
	}
}

} // namespace dawg

#endif //DAWG_PILE_H
