#pragma once
#ifndef DAWG_VECIO_H
#define DAWG_VECIO_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <vector>
#include <iostream>

#include <boost/config.hpp>

namespace dawg {

namespace details {

// vecio_info adapted from boost/tuple/tuple_io.hpp
// Copyright (C) 2001 Jaakko Jarvi (jaakko.jarvi@cs.utu.fi)
//               2001 Gary Powell (gary.powell@sierra.com)
class vecio_info {
public:   
	enum manipulator_type { open, close, delimiter, number_of_manipulators};
private:
	static int get_stream_index (int m) {
		static const int stream_index[number_of_manipulators]
			= { std::ios::xalloc(), std::ios::xalloc(), std::ios::xalloc() };
		return stream_index[m];
	}

	vecio_info(const vecio_info&);
	vecio_info();   

public:
	template<class CharType, class CharTrait>
	static CharType get_manipulator(std::basic_ios<CharType, CharTrait>& i, 
		                   manipulator_type m) {
		CharType c = static_cast<CharType>(i.iword(get_stream_index(m)) ); 
		// parentheses and space are the default manipulators
		if (!c) {
			switch(m) {
			case details::vecio_info::open :  c = i.widen('('); break;
			case details::vecio_info::close : c = i.widen(')'); break;
			case details::vecio_info::delimiter : c = i.widen(' '); break;
			default: break;
			}
		}
		return c;
	}

	template<class CharType, class CharTrait>
	static void set_manipulator(std::basic_ios<CharType, CharTrait>& i, 
		               manipulator_type m, CharType c) {
		i.iword(get_stream_index(m)) = static_cast<long>(c);
	}
};

} //namespace details

template<class CharType>    
class vecio_manipulator {
	const details::vecio_info::manipulator_type mt;
	CharType f_c;
public:
	explicit vecio_manipulator(details::vecio_info::manipulator_type m, 
                             const char c = 0)
     : mt(m), f_c(c) {}
  
template<class CharTrait>
  void set(std::basic_ios<CharType, CharTrait> &io) const {
     details::vecio_info::set_manipulator(io, mt, f_c);
  }
};

template<class CharType>
inline vecio_manipulator<CharType> set_open(const CharType c) {
   return vecio_manipulator<CharType>(details::vecio_info::open, c);
}

template<class CharType>
inline vecio_manipulator<CharType> set_close(const CharType c) {
   return vecio_manipulator<CharType>(details::vecio_info::close, c);
}

template<class CharType>
inline vecio_manipulator<CharType> set_delimiter(const CharType c) {
   return vecio_manipulator<CharType>(details::vecio_info::delimiter, c);
}

template<class CharType, class CharTrait>
inline std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const vecio_manipulator<CharType>& m) {
  m.set(o);
  return o;
}

template<class CharType, class CharTrait, class T, class A>
inline std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const std::vector<T,A> &v) {
	if(!o.good()) return o;
	
	const CharType l = 
	    details::vecio_info::get_manipulator(o, details::vecio_info::open);
	const CharType d = 
		details::vecio_info::get_manipulator(o, details::vecio_info::delimiter);
	const CharType r = 
		details::vecio_info::get_manipulator(o, details::vecio_info::close);
	if(l != 127) o << l;
	typename std::vector<T,A>::const_iterator it = v.begin();
	if(it != v.end())
		o << *(it++);
	if(d != 127)
		for(;it != v.end();++it)
			o << d << *it;
	else
		for(;it != v.end();++it)
			o << *it;
			
	if(r != 127) o << r;
	
	return o;
}

} //namespace dawg
#endif

