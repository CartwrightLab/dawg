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

template<class _T, class _W>
void printnode(typename FingerTree<_T,_W>::node_ptr p) {
	if(p == NULL)
		return;
	cout << "(";
	printnode<_T,_W>(p->left);
	cout << p->val;
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
	
	typedef FingerTree<int,int> FT;
	FT tree;
	FT::node_ptr r = tree.end_node();
	tree.insert(r,1);
	printnode<int,int>(tree.root_node());
	cout << endl;
	tree.insert(r, 2);
	printnode<int,int>(tree.root_node());
	cout << endl;
	tree.insert(r, 3);
	printnode<int,int>(tree.root_node());
	cout << endl;
	tree.insert(r, 4);
	printnode<int,int>(tree.root_node());
	cout << endl;
	tree.insert(r, 5);
	printnode<int,int>(tree.root_node());
	cout << endl;
	tree.insert(r, 6);
	printnode<int,int>(tree.root_node());
	cout << endl;
	tree.insert(r, 7);
	printnode<int,int>(tree.root_node());
	cout << endl;
	tree.insert(r, 8);
	printnode<int,int>(tree.root_node());
	cout << endl;	
	tree.insert(r, 9);
	printnode<int,int>(tree.root_node());
	cout << endl;	
	tree.insert(r, 10);
	printnode<int,int>(tree.root_node());	
	cout << endl;	
	tree.insert(r, 11);
	printnode<int,int>(tree.root_node());
	cout << endl;	
	tree.insert(r, 12);
	printnode<int,int>(tree.root_node());
	cout << endl;	
	tree.insert(r, 13);
	printnode<int,int>(tree.root_node());
	cout << endl;
				
	return 0;
}

ResidueFactory Sequence2::makeSeq;

