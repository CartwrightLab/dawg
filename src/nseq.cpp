#include "nseq.h"

#include <string>
#include <iostream>

ResidueFactory makeSeq;

using namespace std;

int main(int argc, char* argv[]) {

	string ss("ACGTUacgtu");
	double dd[] = {
		1.0, 2.0, 3.0, 4.0, 5.0,
		1.0, 2.0, 3.0, 4.0, 5.0,
	};
	
	ResidueFactory::seq_iterator it = makeSeq(ss.begin(),ss.end(),dd);
	
	ResidueFactory::seq_iterator nit = it;
	for(int i=0;i<ss.size();++i,++nit)
		cout << (int)((nit)->base) << "-" << ((nit)->rate) << "    ";
	cout << endl;
	return 0;
}
