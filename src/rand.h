// rand.h - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_RAND_H
#define DAWG_RAND_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <vector>

typedef boost::mt19937 DawgRng;

extern DawgRng g_rng;
extern boost::uniform_01<DawgRng, double> g_randReal01;

template<class IntType = int, class RealType = double>
class discrete_distribution
{
public:
	typedef IntType result_type;
	typedef RealType input_type;
	typedef std::vector<input_type> storage_type;
	typedef typename storage_type::const_iterator const_iterator;
	
	discrete_distribution() : _P { }

	template<class It>
	discrete_distribution(const It& first, const It& last) : _P(first, last) { }
	
	void reset() { }

	const_iterator begin() const { return _P.begin(); }
	const_iterator end() const { return _P.end(); }

	template<class Engine>
	result_type operator()(Engine& eng)
	{
		input_type r = eng();
		result_type n = result_type(0);
		for(const_iterator cit = _P.begin();
			cit != _P.end() && r < *cit; ++cit)
			++n;
		return n;
	}
	
private:
	storage_type _P;

};

template<class IntType = int, class RealType = double>
class zipf_distribution
{
public:
	typedef IntType result_type;
	typedef RealType input_type;
	explicit zipf_distribution(const input_type& alpha = input_type(3)) : _alpha(alpha)
	{
		assert(alpha > input_type(1));
		init();
	}
	
	void reset() { }

	// Draw from Zipf distribution, with parameter a > 1.0
	// Devroye Luc (1986) Non-uniform random variate generation.
	//     Springer-Verlag: Berlin. p551
	template<class Engine>
	result_type operator()(Engine& eng)
	{
		input_type x,t;
		const input_type one(1);
		do {
		 x = floor(pow(one-eng(), -one/(_alpha-one)));
		 t = pow(one+one/x, _alpha-one);
		} while(eng()*x*(t-one)*_b > t*(_b-one));
		return static_cast<result_type>(x);
	}
	input_type alpha() const { return _alpha; }

protected:
	void init()
	{
		_b = pow(input_type(2), _alpha-input_type(1));
	}
private:
	input_type _alpha;
	input_type _b;
};

template<class IntType = int, class RealType = double>
class truncated_zipf_distribution : public zipf_distribution<IntType, RealType>
{
public:
	typedef IntType result_type;
	typedef RealType input_type;
	typedef zipf_distribution<IntType, RealType> base_type;
	truncated_zipf_distribution(const input_type& alpha = input_type(3),
		const result_type& m = result_type(-1) ) : base_type(alpha), _max(m)
	{
		assert(alpha > input_type(1) && (m == result_type(-1) || m >= result_type(1)));
		init();
	}
	
	void reset() { base_type::reset(); }

	template<class Engine>
	result_type operator()(Engine& eng)
	{
		result_type u = base_type::operator()(eng);
		if(_max == result_type(-1))
			return u;
		while(u > _max)
			u = base_type::operator()(eng);
		return u;
	}

protected:
	void init()
	{
		base_type::init();
	}
private:
	result_type _max;
};


#endif
