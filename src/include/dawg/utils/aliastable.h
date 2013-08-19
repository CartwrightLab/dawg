#pragma once
#ifndef ALIASTABLE_H
#define ALIASTABLE_H

#include <vector>
#include <numeric>
#include <cassert>
#include <ostream>
#include <boost/range.hpp>
#include <boost/cstdint.hpp>

class alias_table {
public:
	typedef boost::uint64_t uint64;
	typedef boost::uint32_t uint32;
	typedef uint32 category_type;
	
	alias_table() { }
	
	template< typename T >
	explicit alias_table(const T &v) {
		create(v);
	}
	
	category_type get(uint64 u) const {
		uint32 x = static_cast<uint32>(u >> shr_);
		uint32 y = static_cast<uint32>(u);
		return ( y < p_[x]) ? x : a_[x];
	}
	
	const std::vector<uint32>& a() const { return a_;}
	const std::vector<uint32>& p() const { return p_;}

	category_type operator()(uint64 u) const {
		return get(u);
	}
	
	// create the alias table
	template< typename T >
	void create(const T &v) {
		std::vector<double> vv(boost::begin(v),boost::end(v));
		create_inplace(vv);
	}

	template< typename T >
	void create(T first, T last) {
		std::vector<double> vv(first,last);
		create_inplace(vv);
	}
	
	// create the alias table
	void create_inplace(std::vector<double> &v) {
		assert(v.size() <= std::numeric_limits<uint32>::max());
		// round the size of vector up to the nearest power of two
		std::pair<std::vector<double>::size_type,int> ru = round_up(v.size());
		const std::vector<double>::size_type sz = ru.first;
		v.resize(sz,0.0);
		a_.resize(sz,0);
		p_.resize(sz,0);
		// use the number of bits to calculate the right shift operand
		shr_ = 64 - ru.second;
		
		// find scale for input vector
		double d = std::accumulate(v.begin(),v.end(),0.0)/sz;
		
		// find first large and small values
		//     g: current large value index
		//     m: current small value index
		//    mm: next possible small value index
		std::vector<double>::size_type g,m,mm;
		for(g=0; g<sz && v[g] <  d; ++g)
			/*noop*/;
		for(m=0; m<sz && v[m] >= d; ++m)
			/*noop*/;
		mm = m+1;
		
		// contruct table
		while(g < sz && m < sz) {
			assert(v[m] < d);
			p_[m] = static_cast<uint32>(4294967296.0/d*v[m]);
			a_[m] = static_cast<uint32>(g);
			v[g] = (v[g]+v[m])-d;
			if(v[g] >= d || mm <= g) {
				for(m=mm; m<sz && v[m] >= d; ++m)
					/*noop*/;
				mm = m+1;
			} else
				m = g;
			for(; g<sz && v[g] <  d; ++g)
				/*noop*/;
		}
		// if we stopped early fill in the rest
		if(g < sz) {
			p_[g] = std::numeric_limits<uint32>::max();
			a_[g] = static_cast<uint32>(g);
			for(g=g+1; g<sz; ++g) {
				if(v[g] < d)
					continue;
				p_[g] = std::numeric_limits<uint32>::max();
				a_[g] = static_cast<uint32>(g);
			}
		}
		// if we stopped early fill in the rest
		if(m < sz) {
			p_[m] = std::numeric_limits<uint32>::max();
			a_[m] = static_cast<uint32>(m);
			for(m=mm; m<sz; ++m) {
				if(v[m] > d)
					continue;
				p_[m] = std::numeric_limits<uint32>::max();
				a_[m] = static_cast<uint32>(m);
			}
		}
	}

	template<class CharType, class CharTrait>
	void print_table(std::basic_ostream<CharType, CharTrait>& o) {
		for(std::size_t n = 0; n < a_.size(); ++n) {
			o << n << "\t" << a_[n] << "\t" << p_[n] << "\n";
		}
	}
	
private:
	template<typename T>
	inline static std::pair<T,int> round_up(T x) {
		T y = static_cast<T>(2);
		int k = 1;
		for(;y < x;y*=2,++k)
			/*noop*/;
		return std::make_pair(y,k);
	}

	std::vector<uint32> a_,p_;
	int shr_;
};

template<class CharType, class CharTrait, class T, class A>
inline std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const alias_table &a) {
	
}


#endif //ALIASTABLE_H
