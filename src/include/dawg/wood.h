#pragma once
#ifndef DAWG_WOOD_H
#define DAWG_WOOD_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace dawg {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

struct wood_node {
	unsigned short anc;
	unsigned short right;
	float length;
	std::string label;
	wood_node() : label(), length(), anc(), right() { }
	wood_node(short r) : label(), length(), anc(), right(r) { }
	wood_node(const std::string &lab, float len=0.0f) : label(lab), length(len), anc(0), right(0) { }
	wood_node(const std::string &lab, const boost::optional<float> &len) :
		label(lab), length(len.get_value_or(0.0f)), anc(0), right(0) {
	}
};
typedef std::vector<wood_node> wood_data;

struct make_inode_impl {
	template<typename V, typename C>
	struct result {typedef void type; };
	template<typename V, typename C>
	void operator()(V &vec, const C& width) const {
		vec.back().anc = 1;
		C w = width+1;
		(vec.end()-w)->anc = w;
		vec.push_back(wood_node(w));
	}
};
const phoenix::function<make_inode_impl> make_inode;

template <typename Iterator>
struct newick_grammar : qi::grammar<Iterator, wood_data&(), ascii::space_type> {
	// http://evolution.genetics.washington.edu/phylip/newick_doc.html
	typedef wood_data::size_type size_type;
	newick_grammar() : newick_grammar::base_type(start) {
		using ascii::space; using boost::spirit::eps;
		using boost::spirit::float_; using boost::spirit::char_;
		using boost::spirit::lexeme; using boost::spirit::raw;
		using boost::spirit::arg_names::_1; using boost::spirit::arg_names::_2;
		using boost::spirit::arg_names::_val; using boost::spirit::arg_names::_r1;
		using phoenix::construct; using phoenix::bind;
		using phoenix::back; using phoenix::push_back;
		
		start    = node(_val) >> ';';
		node     %= tip(_r1) | inode(_r1);
		tip      = label[push_back(_r1,construct<wood_node>(_1))][_val=1]
		           >> -(':' >> float_[bind(&wood_node::length, back(_r1)) = _1]);
		inode    = '(' >> node(_r1)[_val=_1] >> (+(','
		               >> node(_r1)[make_inode(_r1, _1)][_val+=_1+1])
					   | eps[make_inode(_r1,0)][_val+=1]) >> ')'
					>> -(label[bind(&wood_node::label, back(_r1)) = _1] ||
					(':' >> float_[bind(&wood_node::length, back(_r1)) = _1]));
		label    %= unquoted | quoted;
		unquoted %= lexeme[+(char_ - (char_(":,)(;'[]")|space))];
		quoted   %= raw[lexeme['\'' >>
			*(char_ - '\'') >> *(char_("\'") >> char_("\'") >> *(char_ - '\''))
			>> '\'']];
	}
	
	qi::rule<Iterator, wood_data&(), ascii::space_type> start;
	qi::rule<Iterator, size_type(wood_data&), ascii::space_type> node;
	qi::rule<Iterator, size_type(wood_data&), ascii::space_type> tip;
	qi::rule<Iterator, size_type(wood_data&), ascii::space_type> inode;
	qi::rule<Iterator, std::string(), ascii::space_type> label;
	qi::rule<Iterator, std::string(), ascii::space_type> unquoted;
	qi::rule<Iterator, std::string(), ascii::space_type> quoted;	
};

class wood {
public:
	typedef wood_node node;
	wood_data data;
	
	template<typename Iterator>
	bool parse(Iterator first, Iterator last) {
		using ascii::space;
		newick_grammar<Iterator> newick_parser;
		bool r = qi::phrase_parse(first, last, newick_parser, data, space);
		if( first != last || !r )
			return false;
		// The parser produces a tree that is tips first.
		// Transform to root first.
		std::reverse(data.begin(), data.end());
		return true;
	}
};

}

#endif //DAWG_WOOD_H
