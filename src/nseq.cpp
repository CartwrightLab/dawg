#include <iostream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>


#include "sequence.h"

using namespace std;

template<class _T, class _W>
void printnode(typename dawg::finger_tree<_T,_W>::node::pointer p) {
	if(p == NULL)
		return;
	cout << "(";
	printnode<_T,_W>(p->left);
	cout << p->val.base() << ((p->color) ? "r" : "b");
	printnode<_T,_W>(p->right);
	cout << ")";
}

typedef dawg::finger_tree<dawg::residue, dawg::evo_node_weight<> > FT;

int main(int argc, char* argv[]) {
	FT tree;
	char in;
	dawg::residue_factory make_seq;
	
	std::string ss = "ACGTACGACGTACGACGTACG";
	
	make_seq(ss.begin(),ss.end(),tree);
	
	const FT& ctree = tree;
	
	FT::iterator it = tree.find(2ul);
	it->val.base(1);
	
	FT::const_iterator cit = ctree.find(2ul);
	cit->val.base();
	
	return 0;
}

