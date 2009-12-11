#pragma once
#ifndef DAWG_DETAILS_CUMTREE_H
#define DAWG_DETAILS_CUMTREE_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <vector>
#include <limits>

namespace dawg { namespace details {

template<typename _Tp, typename _Alloc = std::allocator<_Tp> >
class cum_tree {
public:
	typedef typename std::vector<_Tp, _Alloc> storage_type;
	typedef typename storage_type::size_type size_type;
	typedef typename storage_type::value_type value_type;

	cum_tree() : mid(0) { }

	cum_tree(size_type n) : data(upper_bound(n), value_type(0)) {
		mid = data.size()/2;
	}

	template<class It>
	cum_tree(It first, It last) {
		assign(first,last);
	}
	
	template<class It>
	void assign(It first, It last) {
		data.assign(upper_bound(last-first),value_type(0));
		mid = data.size()/2;
		typename storage_type::iterator it = data.begin();
		value_type v = value_type(0);
		for(;b != e; ++b) {
			v += *b;
			*(it++) = v;
		}
	}
	
	void increase(size_type u, const value_type& x) {
		for(; u < data.size(); ++u)
			data[u] += x;
	}

	// TODO: Try std::lower_bound
	size_type operator()(const value_type& v) const {
		size_type r = 0;
		for(size_type u = mid; u > 0; u /= 2) {
			if(v >= data[r+u-1])
				r += u;
		}
		return r;
	}

	template<size_type n>
	size_type operator()(const value_type& v) const {
		size_type r = 0;
		for(size_type u = n; u > 0; u /= 2) {
			if(v >= data[r+u-1])
				r += u;
		}
		return r;
	}

protected:
	size_type upper_bound(size_type n) {
		size_type u = 4;
		for(;u<n;u*=2)
			/*noop*/;
		return u;
	}

	storage_type data;
	size_type mid;
};

}}

#endif // DAWG_DETAILS_CUMTREE_H
