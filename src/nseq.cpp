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
	cout << p->val.val << ((p->color) ? "r" : "b");
	printnode<_T,_W>(p->right);
	cout << ")";
}

long unsigned int ul_find = 1ul;
double d_find = 1.1;
dawg::evo_node_weight<> w_find(1.1,1ul);
typedef dawg::finger_tree<dawg::evo_node, dawg::evo_node_weight<> > FT;
extern FT::iterator tit;
FT::iterator tit;

int main(int argc, char* argv[]) {
	FT tree;
	char in;

	for(int i=0;i<256;++i) {
		tree.insert(tree.end(), dawg::evo_node(i));
		/*
		printnode<dawg::evo_node,dawg::evo_node_weight<> >(&*tree.root());
		cout << endl;
		for(FT::iterator it = tree.begin(); it!=tree.end();++it) {
			cout << it->val.val << "/" << it->weight.rate << "/" << it->weight.length << " ";
		}
		cout << endl;
		*/
	}

	for(int i=0;i<1000000000;++i) {
		//tit = tree.find(ul_find);
		tit = tree.find(d_find);
		//tit = tree.find(w_find);
	}		
	return 0;
}

ResidueFactory Sequence2::makeSeq;

