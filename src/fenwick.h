#pragma once
#ifndef DAWG_FENWICK_H
#define DAWG_FENWICK_H

#include <vector>

// Using an alternative implementation of a Fenwick Tree
// http://www.algorithmist.com/index.php/Fenwick_tree
// data[i] = sum( i & (i+1) . . . i )

template<typename _Tp, typename _Alloc = std::allocator<_Tp> >
class fenwick_tree {
public:
	typedef typename std::vector<_Tp, _Alloc> storage_type;
	typedef typename storage_type::size_type size_type;
	typedef typename storage_type::value_type value_type;

	void increase(size_type u, const value_type& x) {
		if(u >= data.size())
			data.resize(u+1);
		for(; u < data.size(); u |= u+1)
			data[u] += x;
	}
	size_type operator()(const value_type& v) {


	}


protected:
	void resize(size_type x) {
		size_type z = data.size();
		if(z == 0)
			z = 1;
		while(z < x)
			z << 1;
		data.resize(z, value_type(0));
		mask = (z-1) >> 1;
	}

	storage_type data;
	size_type mask;
};

#endif // DAWG_FENWICK_H
