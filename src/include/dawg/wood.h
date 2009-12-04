#pragma once
#ifndef DAWG_WOOD_H
#define DAWG_WOOD_H

/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/scoped_ptr.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_variable.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/void.hpp>
#include <boost/variant/recursive_variant.hpp>

namespace dawg {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

//typedef std::pair<std::string,double> wood_tip;
struct wood_tip {
	wood_tip() : label(), length() {}
	std::string label;
	double length;
	
	/*wood_tip& operator=(const wood_tip& r) {
		label = r.label;
		length = r.length;
		std::cerr << '*';
		return *this;
	}
	*/
};

struct wood_inode;
typedef boost::variant<wood_tip, boost::recursive_wrapper<wood_inode> > wood_node;
struct wood_inode {
	std::vector<wood_node> children;
	wood_tip self;
};

}

BOOST_FUSION_ADAPT_STRUCT( dawg::wood_tip,
	(std::string, label)(double, length)
)
BOOST_FUSION_ADAPT_STRUCT( dawg::wood_inode,
	(std::vector<dawg::wood_node>, children)(dawg::wood_tip, self)
)

namespace dawg{

void tipprint(const wood_tip& t) {
	std::cout << t.label << " " << t.length << std::endl;
}

template <typename Iterator>
struct newick_grammar : qi::grammar<Iterator, dawg::wood_node(), ascii::space_type> {
	// http://evolution.genetics.washington.edu/phylip/newick_doc.html

	newick_grammar() : newick_grammar::base_type(start) {
		using boost::spirit::double_;
		using boost::spirit::lexeme;
		using boost::spirit::char_;
		using ascii::space;
		using boost::spirit::arg_names::_1;
		using boost::spirit::arg_names::_val;
		using boost::phoenix::arg_names::arg1;
		using boost::spirit::raw;
		using phoenix::val;
		using phoenix::for_each;
		using phoenix::bind;
		
		start    %= node >> ';';
		node     %= tip | inode;
		tip      %= label >> -(':' >> double_);
		tip2     %= label || (':' >> double_);
		inode    %= '(' >> (node % ',') >> ')' >> -tip2;
		label    %= unquoted | quoted;
		unquoted %= lexeme[+(char_ - (char_(":,)(;'[]")|space))];
		quoted   %= raw[lexeme['\'' >>
			*(char_ - '\'') >> *(char_("\'") >> char_("\'") >> *(char_ - '\''))
			>> '\'']];
	}
	
	qi::rule<Iterator, dawg::wood_node(), ascii::space_type> start;
	qi::rule<Iterator, dawg::wood_node(), ascii::space_type> node;
	qi::rule<Iterator, dawg::wood_inode(), ascii::space_type> inode;
	qi::rule<Iterator, dawg::wood_tip(), ascii::space_type> tip;
	qi::rule<Iterator, dawg::wood_tip(), ascii::space_type> tip2;
	qi::rule<Iterator, std::string(), ascii::space_type> label;
	qi::rule<Iterator, std::string(), ascii::space_type> unquoted;
	qi::rule<Iterator, std::string(), ascii::space_type> quoted;
};

// Represents a phylogenetic tree structure
class wood {
public:
	wood_node root;
	
	template<typename Iterator>
	bool parse(Iterator first, Iterator last) {
		using ascii::space;
		
		newick_grammar<Iterator> newick_parser;
		
		bool r = qi::phrase_parse(first, last, newick_parser, root, space);
		if( first != last )
			return false;
		return r;
	}
};

}
#endif //DAWG_WOOD_H

