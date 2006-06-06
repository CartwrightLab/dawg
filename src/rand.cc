// rand.cc - Copyright (C) 2004-2006 Reed A. Cartwright (all rights reserved)

#include "rand.h"

DawgRng g_rng;
//boost::uniform_01<DawgRng, double> g_rand01(g_rng);
boost::uniform_real<> uni_dist(0,1);
boost::variate_generator<DawgRng&, boost::uniform_real<> > g_randReal01(g_rng, uni_dist);
zipf_distribution<> zipf_dist(1.5);
boost::variate_generator<DawgRng&, zipf_distribution<> > g_randZipf(g_rng, zipf_dist);

//boost::gamma_distribution<> gamma_dist();
//boost::variate_generator<DawgRng&, boost::gamma_distribution<> > g_randGamma(g_rng, gamma_dist);

// Draw from Zipf distribution, with parameter a > 1.0
// Devroye Luc (1986) Non-uniform random variate generation.
//     Springer-Verlag: Berlin. p551
//unsigned int rand_zipf(double a)
//{
//	double b = pow(2.0, a-1.0);
//	double x,t;
//	do {
//	 x = floor(pow(1.0-g_randReal01(), -1.0/(a-1.0)));
//	 t = pow(1.0+1.0/x, a-1.0);
//	} while( g_randReal01()*x*(t-1.0)*b > t*(b-1.0));
//	return static_cast<unsigned int>(x);
//}

