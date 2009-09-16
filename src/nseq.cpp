#include <sys/time.h>

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

long unsigned int ul_find[256];
double d_find[256];
dawg::evo_node_weight<> w_find[256];
typedef dawg::finger_tree<dawg::residue, dawg::evo_node_weight<> > FT;
//typedef dawg::finger_tree<dawg::evo_node, weight> FT;
extern FT::iterator tit;
FT::iterator tit;

int main(int argc, char* argv[]) {
	FT tree;
	char in;
	
	for(int i=0;i<50;++i) {
		tree.insert(tree.end(), dawg::residue());
		ul_find[i] = i;
		d_find[i] = i;
		w_find[i] = dawg::evo_node_weight<>(i,i);
	}
		
	for(int i=0;i<50;++i) {
		tit = tree.find(ul_find[i]);
		cout << tit->val.base() << " ";
		tit = tree.find(d_find[i]);
		cout << tit->val.base() << " ";
		tit = tree.find(w_find[i]);
		cout << tit->val.base() << endl;
	}
		
	return 0;
}

ResidueFactory Sequence2::makeSeq;

