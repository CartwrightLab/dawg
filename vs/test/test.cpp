// test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "rand.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	zipf_distribution<> zipf_dist(1.5);
	boost::variate_generator<DawgRng&, zipf_distribution<> > rand_zipf(g_rng, zipf_dist);
	for(int i=0;i<300;++i)
		cout << rand_zipf() << " ";
	cout << endl;
	return 0;
}

