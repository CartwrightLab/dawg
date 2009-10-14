#pragma once
#ifndef DAWG_BITREE_H
#define DAWG_BITREE_H

#include <vector>
#include <limits>

// Using an alternative implementation of a binary indexed tree
// http://www.algorithmist.com/index.php/Fenwick_tree
// data[i] = sum( 0 . . . i )

template<typename _Tp, typename _Alloc = std::allocator<_Tp> >
class bitree {
public:
	typedef typename std::vector<_Tp, _Alloc> storage_type;
	typedef typename storage_type::size_type size_type;
	typedef typename storage_type::value_type value_type;

	bitree() : mid(0) { }

	bitree(size_type n) : data(upper_bound(n), value_type(0)) {
		mid = data.size()/2;
	}

	template<class It>
	bitree(It b, It e) : data(upper_bound(e-b), value_type(0)) {
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

#endif // DAWG_FENWICK_H
