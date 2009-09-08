#include <iostream>

#include "nseq.h"

#include <string>

ResidueFactory makeSeq;

using namespace std;

void printseq(const Sequence2 &seq2) {
	cout << seq2.size() << ": ";
	
	for(Sequence2::const_iterator nit = seq2.begin(); nit != seq2.end(); ++nit)
		cout << (int)((nit)->base) << "-" << ((nit)->rate) << "    ";
	cout << endl;
	for(Sequence2::Parts::const_iterator it = seq2._parts().begin(); it != seq2._parts().end(); ++it)
		cout << it->size() << " ";
	cout << endl;
}

int main(int argc, char* argv[]) {
	string ss("ACGT");
	double dd[] = {
		1.0, 2.0, 3.0, 4.0,
	};
	
	Sequence2 seq2;
	seq2.append(ss.begin(),ss.end(),dd);
	printseq(seq2);
	
	Sequence2 seq2a = seq2;
	printseq(seq2a);
	
	seq2.insert(3, seq2a);
	seq2.insert(8, seq2a);
	printseq(seq2);
	seq2.remove(1, 4);
	printseq(seq2);	
	seq2.remove(1, 4);
	printseq(seq2);
	
	return 0;
}

ResidueFactory Sequence2::makeSeq;

