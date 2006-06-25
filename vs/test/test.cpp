// test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "time.h"
#include <iostream>
#include "../../src/sequence.h"
#include "../../src/indelmodel.h"

using namespace std;
using namespace boost;

void print_seq(const Dawg::Sequence& seq)
{
	for(int i = 0; i < seq.Length(); ++i)
		cout << seq.Base(i);
}

int _tmain(int argc, _TCHAR* argv[])
{
	g_rng.seed(time(NULL));
	std::vector<double> freqs(4, 0.25);
	discrete_distribution<unsigned int> base_dist(freqs.begin(), freqs.end());
	Dawg::SequenceFactory fac(base_dist);
	Dawg::Sequence seq = fac(50);
	print_seq(seq);
	cout << endl;

	Dawg::GillespieProcessor gp;
	truncated_zipf_distribution<unsigned int> zd(1.5, 5);
	gp.AddElement(new Dawg::Deletion<truncated_zipf_distribution<unsigned int> >(0.1, zd));
	gp.AddElement(new Dawg::Insertion<truncated_zipf_distribution<unsigned int> >(0.1, zd, fac));

	gp(seq, 0.1);
	print_seq(seq);
	cout << endl;

	return 0;
}

