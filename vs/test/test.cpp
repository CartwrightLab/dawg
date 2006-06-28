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

void print_aln(const Dawg::Sequence& seq)
{
	for(Dawg::Sequence::aln_seq_type::const_iterator cit = seq.Aln().begin();
		cit != seq.Aln().end(); ++cit)
		cout << *cit;
}

template<class T, class A, size_t N>
const T& assign_from_array(T& t, const A (&a)[N])
{
	t.assign(&a[0], &a[N]);
	return t;
}

int _tmain(int argc, _TCHAR* argv[])
{
	vector<double> v;
	double a[] = {1, 2, 3, 4, 5,1, 2, 3, 4, 5,1, 2, 3, 4, 5};
	assign_from_array(v, a);
	for(vector<double>::iterator it = v.begin(); it != v.end(); ++it)
		cout << *it << endl;
	
	//g_rng.seed(time(NULL));
	//std::vector<double> freqs(4, 0.25);
	//discrete_distribution<unsigned int> base_dist(freqs.begin(), freqs.end());
	//Dawg::SequenceFactory fac(base_dist);
	//Dawg::Sequence seq = fac(10);
	//Dawg::Sequence seq2 = fac(10);
	//Dawg::Sequence seq3 = fac(5);
	//seq.AppendSection(seq2);
	//seq.Insert(10, seq3);
	//seq.Insert(25, seq3);
	//seq.Insert(0, seq3);
	//seq.Delete(0, 20);

	//print_seq(seq);
	//cout << endl;
	//print_aln(seq);
	//cout << endl;

	//seq.ReadSection(0, seq3);
	//print_seq(seq3);
	//cout << endl;
	//print_aln(seq3);
	//cout << endl;

	//seq.ReadSection(1, seq3);
	//print_seq(seq3);
	//cout << endl;
	//print_aln(seq3);
	//cout << endl;

	////Dawg::GillespieProcessor gp;
	////truncated_zipf_distribution<unsigned int> zd(1.5, 5);
	////gp.AddElement(new Dawg::Deletion<truncated_zipf_distribution<unsigned int> >(0.1, zd));
	////gp.AddElement(new Dawg::Insertion<truncated_zipf_distribution<unsigned int> >(0.1, zd, fac));

	////gp(seq, 0.25);
	////print_seq(seq);
	////cout << endl << endl;
	////print_aln(seq);
	////cout << endl << endl;

	return 0;
}

