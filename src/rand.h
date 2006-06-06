// rand.h - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_RAND_H
#define DAWG_RAND_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>

template<class IntType = unsigned int, class RealType = double>
class zipf_distribution
{
public:
	typedef IntType result_type;
	typedef RealType input_type;
	explicit zipf_distribution(const input_type& alpha) : _alpha(alpha)
	{
		assert(alpha > input_type(1));
		init();
	}
	
	void reset() { }

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

private:
	void init()
	{
		_b = pow(input_type(2), _alpha-input_type(1));
	}
	boost::uniform_real<RealType> _uni;
	input_type _alpha;
	input_type _b;
};

typedef boost::mt19937 DawgRng;

extern DawgRng g_rng;
//extern boost::uniform_01<DawgRng, double> g_rand01;
extern boost::variate_generator<DawgRng&, boost::uniform_real<> > g_randReal01;
extern boost::variate_generator<DawgRng&, boost::gamma_distribution<> > g_randGamma;
extern boost::variate_generator<DawgRng&, zipf_distribution<> > g_randZipf;

#endif
