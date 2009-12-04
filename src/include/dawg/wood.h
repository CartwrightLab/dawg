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
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/include/phoenix_function.hpp>



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

namespace dawg { namespace ex {

struct wood_node {
	unsigned short anc;
	unsigned short right;
	float length;
	std::string label;
	wood_node() : label(), length(), anc(), right() { }
	wood_node(const std::string &lab, float len=0.0f) : label(lab), length(len), anc(0), right(0) { }
	wood_node(const std::string &lab, const boost::optional<float> &len) :
		label(lab), length(len.get_value_or(0.0f)), anc(0), right(0) {
	}
};

typedef std::vector<wood_node> wood_data;

struct make_tip_impl {
	template<typename V, typename L, typename N>
	struct result {
		typedef void type;
	};
	
	template<typename V, typename L, typename N>
	void operator()(V &val, const L& lab, const N& len) {
		val.assign(1, wood_node(lab,len));
	}
};

phoenix::function<make_tip_impl> make_tip;

struct make_node_impl {
	template<typename V, typename C>
	struct result {typedef void type; };
	template<typename V, typename C>
	void operator()(V &val, const C& cont) const {
		val.assign(1, wood_node());
		val.reserve(256); //yeah we could be smarter about this, but this is simpler.
		typename C::size_type z = cont.size();
		if(z == 2) {
			val.insert(val.end(), cont[0].begin(), cont[0].end());
			val[1].anc = 1;
			val.insert(val.end(), cont[1].begin(), cont[1].end());
			typename C::value_type::size_type r = 1+cont[0].size();
			val[r].anc = r;
			val[0].right = r;
		} else if(z > 2) {
			val.insert(val.end(), cont[0].begin(), cont[0].end());
			val[1].anc = 1;
			typename C::value_type::size_type r = 1+cont[0].size(),p=0;
			
			for(typename C::size_type x = 1; x < z-1; ++x) {
				val.push_back(wood_node());
				val[p].right = r;
				val[p+=r].anc = r;
				val.insert(val.end(), cont[x].begin(), cont[x].end());
				val[p+1].anc = 1;
				r = 1+cont[x].size();
			}
			val.insert(val.end(), cont[z-1].begin(), cont[z-1].end());
			val[p+r].anc = r;
			val[p].right = r;
		} else if(z == 1) {
			val.insert(val.end(), cont[0].begin(), cont[0].end());
			val[0].right = 1;
			val[1].anc = 1;
		}
	}
};
phoenix::function<make_node_impl> make_node;

template <typename Iterator>
struct newick_grammar2 : qi::grammar<Iterator, wood_data(), ascii::space_type> {
	// http://evolution.genetics.washington.edu/phylip/newick_doc.html
	newick_grammar2() : newick_grammar2::base_type(start) {
		using boost::spirit::float_;
		using boost::spirit::lexeme;
		using boost::spirit::char_;
		using ascii::space;
		using boost::spirit::arg_names::_1;
		using boost::spirit::arg_names::_2;
		using boost::spirit::arg_names::_val;
		using boost::phoenix::arg_names::arg1;
		using boost::phoenix::construct;
		using boost::spirit::raw;
		using phoenix::val;
		using phoenix::assign;
		using phoenix::for_each;
		using phoenix::bind;
		
		start    %= node >> ';';
		node     = (label[assign(_val, 1, construct<wood_node>(_1))] >> -(':'
					>> float_[bind(&wood_node::length, _val[0]) = _1]
				 ))
		         | ('(' >> (node % ',')[make_node(_val,_1)] >> ')'
				     >> -(label[bind(&wood_node::label, _val[0]) = _1]
				     || (':' >> float_[bind(&wood_node::length, _val[0]) = _1]
				 )));
		//tip      = (label >> -(':' >> float_))[assign(_val, 1, construct<wood_node>(_1,_2))];
		label    %= unquoted | quoted;
		unquoted %= lexeme[+(char_ - (char_(":,)(;'[]")|space))];
		quoted   %= raw[lexeme['\'' >>
			*(char_ - '\'') >> *(char_("\'") >> char_("\'") >> *(char_ - '\''))
			>> '\'']];
	}
	
	qi::rule<Iterator, wood_data(), ascii::space_type> start;
	qi::rule<Iterator, wood_data(), ascii::space_type> node;
	qi::rule<Iterator, std::string(), ascii::space_type> label;
	qi::rule<Iterator, std::string(), ascii::space_type> unquoted;
	qi::rule<Iterator, std::string(), ascii::space_type> quoted;
};

class wood {
public:
	wood_data data;
	
	template<typename Iterator>
	bool parse(Iterator first, Iterator last) {
		using ascii::space;
		
		newick_grammar2<Iterator> newick_parser;
		
		bool r = qi::phrase_parse(first, last, newick_parser, data, space);
		if( first != last )
			return false;
		return r;
	}
};


}}

#endif //DAWG_WOOD_H

