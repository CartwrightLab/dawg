#pragma once
#ifndef DAWG_WOOD_H
#define DAWG_WOOD_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/spirit/include/version.hpp>
#if SPIRIT_VERSION < 0x2010
#	error Spirit version 2.1 or greater required.
#endif

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <set>
#include <iterator>


namespace dawg {

namespace qi = boost::spirit::qi;
namespace standard = boost::spirit::standard;
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
	inline bool unlabeled() const {
		return (label.empty() || label[0] == '(');
	}
	inline bool terminal() const {
		return (right == 0);
	}
	inline bool one_child() const {
		return (right == 1);
	}
	inline bool two_child() const {
		return (right > 1);
	}
};

class wood {
public:
	typedef wood_node node;
	typedef std::vector<node> data_type;

	inline const data_type& data() const { return _data; }
	inline const std::string& root_label() const { return root_name; }
	inline bool  has_desc(const std::string &d) const  { return desc_names.count(d) != 0; }
	inline const std::set<std::string>& descs() const { return desc_names; }

	template<typename Iterator>
	bool parse(Iterator first, Iterator last);
	
	template<std::size_t _N>
	bool parse(const char (&str)[_N]) {
		return parse(&str[0], &str[_N]);
	}
	bool parse(const std::string &str) {
		return parse(str.begin(), str.end());
	}
	
	bool autolabel() {
		if(_data.empty())
			return true;
		desc_names.clear();
		for(data_type::reverse_iterator it = _data.rbegin(); it != _data.rend(); ++it) {
			if(it->unlabeled() && !it->terminal()) {
				if(it->one_child()) {
					it->label = "(" + get_right(it)->label + ")";
				} else {
					std::string &a = get_left(it)->label;
					std::string &b = get_right(it)->label;
					if(b > a)
						it->label = "(" + a + "," + b + ")";
					else
						it->label = "(" + b + "," + a + ")";
				}
			}
			// add label to descendents list
			if(!desc_names.insert(it->label).second)
				return DAWG_ERROR("invalid tree; node label '" << it->label
					              << "' used more than once by Tree.Tree.");
		}
	
		// set root name and remove it from the desc set
		root_name = _data.front().label;
		desc_names.erase(root_name);
		return true;
	}
	
	
	template<typename Iterator>
	static Iterator get_left(Iterator me) {
		std::advance(me, 1);
		return me;
	}
	template<typename Iterator>
	static Iterator get_right(Iterator me) {
		std::advance(me, static_cast<int>(me->right));
		return me;
	}
	template<typename Iterator>
	static Iterator get_anc(Iterator me) {
		std::advance(me, -static_cast<int>(me->anc));
		return me;
	}
	template<typename Iterator>
	static std::reverse_iterator<Iterator> get_left(std::reverse_iterator<Iterator> me) {
		std::advance(me, -1);
		return me;
	}
	template<typename Iterator>
	static std::reverse_iterator<Iterator> get_right(std::reverse_iterator<Iterator> me) {
		std::advance(me, -static_cast<int>(me->right));
		return me;
	}
	template<typename Iterator>
	static std::reverse_iterator<Iterator> get_anc(std::reverse_iterator<Iterator> me) {
		std::advance(me, static_cast<int>(me->anc));
		return me;
	}	

protected:
	data_type _data;
	std::set<std::string> desc_names;
	std::string root_name;
	
};


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
		           >> -(':' >> float_[bind(&wood_node::length, back(_r1)) = _1]);
		inode    = '(' >> node(_r1)[_val=_1] >> (+(','
		               >> node(_r1)[make_inode(_r1, _1)][_val+=_1+1])
					   | eps[make_inode(_r1,0)][_val+=1]) >> ')'
					>> -(label[bind(&wood_node::label, back(_r1)) = _1] ||
					(':' >> float_[bind(&wood_node::length, back(_r1)) = _1]));
		label    = unquoted | quoted;
		unquoted = lexeme[+(char_ - (char_(":,)(;'[]")|space))];
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

}

#endif //DAWG_WOOD_H
