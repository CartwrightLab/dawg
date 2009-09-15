#include <iostream>
#include <string>

#include "nseq.h"

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

template<class _T, class _W>
void printnode(typename dawg::finger_tree<_T,_W>::node::pointer p) {
	if(p == NULL)
		return;
	cout << "(";
	printnode<_T,_W>(p->left);
	cout << p->val << ((p->color) ? "r" : "b");
	printnode<_T,_W>(p->right);
	cout << ")";
}

int main(int argc, char* argv[]) {
	/*
	string ss("ACGT");
	double dd[] = {
		1.0, 2.0, 3.0, 4.0,
	};
	
	Sequence2 seq2;
	seq2.append(ss.begin(),ss.end(),dd);
	printseq(seq2);
	
	Sequence2 seq2a = seq2.clone();
	printseq(seq2a);
	
	seq2.insert(3, seq2a.clone());
	seq2.insert(8, seq2a.clone());
	printseq(seq2);
	seq2a.begin()->base = 10;
	printseq(seq2);
	seq2.erase(1, 4);
	printseq(seq2);	
	seq2.erase(1, 4);
	printseq(seq2);
	seq2.optimize();
	printseq(seq2);
	*/
	
	typedef dawg::finger_tree<int,int> FT;
	FT tree;
	char in;
	for(int i=21;i<28;++i) {
		tree.insert(tree.begin(), i);
		printnode<int,int>(&*tree.root());
		cout << endl;
		for(FT::iterator it = tree.begin(); it!=tree.end();++it) {
			cout << it->val << " ";
		}
		for(FT::iterator it = tree.end(); it!=tree.begin();)
			cout << (--it)->val << " ";
		cout << endl;			
	}
	for(int i=11;i<19;++i) {
		tree.insert(tree.end(), i);
		printnode<int,int>(&*tree.root());
		cout << endl;
	}

	for(FT::iterator it = tree.begin(); it!=tree.end();++it) {
		cout << it->val << " " << it->up->val << " ";
		cin.get(in);
	}
	cout << endl;

	for(FT::iterator it = tree.end(); it!=tree.begin();)
		cout << (--it)->val << " ";
	cout << endl;

	
	return 0;
}

ResidueFactory Sequence2::makeSeq;

