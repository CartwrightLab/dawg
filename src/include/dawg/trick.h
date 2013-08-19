#pragma once
#ifndef DAWG_TRICK_H
#define DAWG_TRICK_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <utility>
#include <string>
#include <vector>
#include <map>

#include <dawg/log.h>
#include <dawg/utils/foreach.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/erase.hpp>

namespace dawg {

struct trick {
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
		
		//inline void get(const std::string& k, trick::section& r) const;
		
		inline void read_aliases();
		
	private:
		static inline void conv(const std::string& ss, std::string& r);
		static inline void conv(const std::string& ss, double& r);
		static inline void conv(const std::string& ss, bool& r);
		static inline void conv(const std::string& ss, unsigned int& r);
		static inline void conv(const std::string& ss, int& r);
		
		inline void read_alias(const std::string& a, const std::string& b);
	};
	typedef std::vector<section> data_type;
	data_type data;

	static bool parse_file(trick &p, const char *cs);
	template<typename Iterator>
	bool parse(Iterator first, Iterator last);
	template<typename Char, typename Traits>
	inline bool parse_stream(std::basic_istream<Char, Traits>& is);
	
	trick() {
		data.push_back(section());
		data.back().name = "_initial_";
		data.back().inherits = "_default_";
	}
	
	inline void read_aliases();	
};

template<typename T>
inline void trick::section::get(const std::string& k, T& r) const {
	db_type::const_iterator it;
	if((it = db.find(k)) != db.end() && !it->second.empty()) {
		section::conv(it->second.front(), r);
	}
}

template<typename T, typename A>
inline void trick::section::get(const std::string& k, std::vector<T,A>& r) const {
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

template<>
inline void trick::section::get(const std::string& k, trick::section& r) const {
	using boost::algorithm::starts_with;
	using boost::algorithm::erase_head_copy;
	r.name = k;
	r.inherits = "_nothing_";
	r.db.clear();
	db_type::const_iterator first = db.lower_bound(k);
	db_type::const_iterator last;
	for(last = first; last != db.end()
		&& starts_with(last->first, k); ++last) {
		r.db.insert(r.db.end(),
			make_pair(erase_head_copy(last->first, static_cast<int>(k.length())), last->second));
	}
}


inline void trick::section::conv(const std::string& ss, std::string& r) {
	r = ss;
}

inline void trick::section::conv(const std::string& ss, double& r) {
	r = strtod(ss.c_str(), NULL);
}

inline void trick::section::conv(const std::string& ss, unsigned int& r) {
	r = strtoul(ss.c_str(), NULL, 0);
}

inline void trick::section::conv(const std::string& ss, int& r) {
	r = strtol(ss.c_str(), NULL, 0);
}

// A value is false if it is equal to 0, f, false, off, no, or blank
inline void trick::section::conv(const std::string& ss, bool& r) {
	using boost::algorithm::iequals;
	r = !(ss.empty() || iequals(ss, "false") || iequals(ss, "0") || iequals(ss, "f")
		|| iequals(ss, "off") || iequals(ss, "no"));
}

inline void trick::section::read_alias(const std::string& a, const std::string& b) {
	db_type::const_iterator it;
	// if b exists or a doesn't, stop
	if(db.find(b) != db.end() || (it = db.find(a)) == db.end())
		return;
	db[b] = it->second;
}

inline void trick::section::read_aliases() {
	#define XM(aname, bname) read_alias(XP(aname), XP(bname));
	#include <dawg/details/aliases.xmh>
	#undef XM
}

inline void trick::read_aliases() {
	foreach(section &sec, data) {
		sec.read_aliases();
	}
}

} // namespace dawg

#endif //DAWG_TRICK_H
