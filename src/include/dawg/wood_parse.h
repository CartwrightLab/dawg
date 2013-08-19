#pragma once
#ifndef DAWG_WOOD_PARSE_H
#define DAWG_WOOD_PARSE_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/object/construct.hpp>
#include <boost/phoenix/bind/bind_member_variable.hpp>
#include <boost/phoenix/stl/container.hpp>
#include <boost/phoenix/operator.hpp>

#include <boost/spirit/include/version.hpp>
#if SPIRIT_VERSION < 0x2020
#	error Spirit version 2.2 or greater required.
#endif

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable: 4127 )
#endif

#include <boost/spirit/include/qi.hpp>

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include<iterator>

#include<dawg/wood.h>

namespace dawg {


namespace qi = boost::spirit::qi;
namespace standard = boost::spirit::standard;
namespace phoenix = boost::phoenix;

struct make_inode_impl {
	typedef void result_type;

	template<typename V, typename C>
	void operator()(V &vec, const C& width) const {
		vec.back().anc = 1;
		wood_node::id_t w = static_cast<wood_node::id_t>(width+1);
		(vec.end()-w)->anc = w;
		vec.push_back(wood_node(w));
	}
};
const phoenix::function<make_inode_impl> make_inode;

template <typename Iterator>
struct newick_grammar : qi::grammar<Iterator, wood::data_type(), standard::space_type> {
	// http://evolution.genetics.washington.edu/phylip/newick_doc.html
	typedef wood::data_type::size_type size_type;
	newick_grammar() : newick_grammar::base_type(start) {
		using standard::space; using standard::char_; using qi::eps;
		using qi::float_; using qi::lexeme; using qi::raw; using qi::omit;
		using qi::_1; using qi::_2; using qi::_val; using qi::_r1;
		using phoenix::construct; using phoenix::bind;
		using phoenix::back; using phoenix::push_back;
		using phoenix::size;
		
		start    = omit[node(_val)] >> ';';
		node     = tip(_r1) | inode(_r1);
		tip      = label[push_back(_r1,construct<wood_node>(_1))][_val=1]
		           >> -(':' >> float_[phoenix::bind(&wood_node::length, back(_r1)) = _1]);
		inode    = '(' >> node(_r1)[_val=_1] >> (+(','
		               >> node(_r1)[make_inode(_r1, _1)][_val+=_1+1])
					   | eps[make_inode(_r1,0)][_val+=1]) >> ')'
					   >> -(label[phoenix::bind(&wood_node::label, back(_r1)) = _1] ||
					   (':' >> float_[phoenix::bind(&wood_node::length, back(_r1)) = _1]));
		label    = unquoted | quoted;
		// Due to the way hidden nodes are constructed,
		// unquoted labels should not begin with {, |, or }.
		unquoted = lexeme[(char_ - (char_(":,)(;'[]|{}")|space)) >>
		                 *(char_ - (char_(":,)(;'[]")|space))];
		quoted   = raw[lexeme['\'' >>
			*(char_ - '\'') >> *(standard::string("\'\'") >> *(char_ - '\''))
			>> '\'']];
	}
	
	qi::rule<Iterator, wood::data_type(), standard::space_type> start;
	qi::rule<Iterator, size_type(wood::data_type&), standard::space_type> node;
	qi::rule<Iterator, size_type(wood::data_type&), standard::space_type> tip;
	qi::rule<Iterator, size_type(wood::data_type&), standard::space_type> inode;
	qi::rule<Iterator, std::string(), standard::space_type> label;
	qi::rule<Iterator, std::string(), standard::space_type> unquoted;
	qi::rule<Iterator, std::string(), standard::space_type> quoted;	
};

template<typename Iterator>
bool wood::parse(Iterator first, Iterator last) {
	using standard::space;
	newick_grammar<Iterator> newick_parser;
	bool r = qi::phrase_parse(first, last, newick_parser, space, _data);
	if( first != last || !r )
		return false;
				
	// The parser produces a tree that is tips first.
	// Transform to root first.
	std::reverse(_data.begin(), _data.end());
	
	if(!autolabel())
		return false;
	
	return true;
}

} // namespace dawg
#endif //DAWG_WOOD_PARSE_H
