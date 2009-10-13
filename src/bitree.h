#pragma once
#ifndef DAWG_BITREE_H
#define DAWG_BITREE_H

#include <vector>
#include <limits>

// Using an alternative implementation of a binary indexed tree
// http://www.algorithmist.com/index.php/Fenwick_tree
// data[i] = sum( i & (i+1) . . . i )

template<typename _Tp, typename _Alloc = std::allocator<_Tp> >
class bitree {
public:
	typedef typename std::vector<_Tp, _Alloc> storage_type;
	typedef typename storage_type::size_type size_type;
	typedef typename storage_type::value_type value_type;

	bitree(size_type n) : data(n, value_type(0)) {
		mid = mk_mid(n);
	}


	// b = (b & (b + 1)) - 1
	void increase(size_type u, const value_type& x) {
		for(; u < data.size(); u |= u+1)
			data[u] += x;
	}
	size_type operator()(const value_type& v) {
		size_type r = 0, u = mask;
		for(size_type u = mask; u > 0 && r < data.size(); u >>= 1) {
			if(v < data[r+u])
				continue;
			v -= data[r+u];
			r += u+1;
		}
		if(v >= data[r])
			++r;
		return (r < data.size() ? r : data.size());
	}

protected:
	size_type mk_mid(size_type n) {
		n = (n-1)/2;
		for(int x=1; x < std::numeric_limits<size_type>::digits; x <<= 1)
			n |= (n >> x);
		return n;
	}

	storage_type data;
	size_type mid;
};

#endif // DAWG_FENWICK_H
