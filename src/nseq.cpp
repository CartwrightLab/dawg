#include <iostream>

#include "nseq.h"

#include <string>

ResidueFactory makeSeq;

using namespace std;

int main(int argc, char* argv[]) {

	string ss("ACGTUacgtu");
	double dd[] = {
		1.0, 2.0, 3.0, 4.0, 5.0,
		1.0, 2.0, 3.0, 4.0, 5.0,
	};
	
	Sequence2 seq2;
	seq2.append(ss.begin(),ss.end(),dd);
	
	

	Sequence2::iterator nit;
	for(nit = seq2.begin(); nit != seq2.end(); ++nit)
		cout << (int)((nit)->base) << "-" << ((nit)->rate) << "    ";
	cout << endl;
	return 0;
}

ResidueFactory Sequence2::makeSeq;
