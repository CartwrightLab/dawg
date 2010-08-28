#pragma once
#ifndef DAWG_WOOD_H
#define DAWG_WOOD_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <set>
#include <iterator>
#include <functional>

#include <boost/optional.hpp>

namespace dawg {

struct wood_node {
	std::string label;
	float length;
	unsigned short anc;
	unsigned short right;

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
	inline const std::set<std::string>& desc_labels() const { return desc_names; }

	template<typename Iterator>
	bool parse(Iterator first, Iterator last);
	
	template<std::size_t _N>
	bool parse(const char (&str)[_N]) {
		return parse(&str[0], &str[_N]);
	}
	bool parse(const std::string &str) {
		return parse(str.begin(), str.end());
	}
	static bool parse_string(wood &w, const std::string &str);

	inline void scale(double d) {
		for(data_type::iterator it = _data.begin(); it != _data.end(); ++it) {
			it->length *= static_cast<float>(d);
		}
	}
	
	bool autolabel() {
		if(_data.empty())
			return true;
		desc_names.clear();
		for(data_type::reverse_iterator it = _data.rbegin(); it != _data.rend(); ++it) {
			if(it->unlabeled() && !it->terminal()) {
				if(it->one_child()) {
					it->label = "{" + get_right(it)->label + "}";
				} else {
					std::string &a = get_left(it)->label;
					std::string &b = get_right(it)->label;
					if(b > a)
						it->label = "{" + a + "," + b + "}";
					else
						it->label = "{" + b + "," + a + "}";
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

}

#endif //DAWG_WOOD_H
