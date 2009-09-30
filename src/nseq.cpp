#include <iostream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>


#include "nseq.h"

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

long unsigned int ul_find[256];
double d_find[256];
dawg::evo_node_weight<> w_find[256];
typedef dawg::finger_tree<dawg::residue, dawg::evo_node_weight<> > FT;
extern FT::iterator tit;
FT::iterator tit;
extern ptime tstart,tend;
ptime tstart,tend;

vector<dawg::residue> store;

int main(int argc, char* argv[]) {
	FT tree;
	char in;
	dawg::residue_factory make_seq;
	
	std::string ss = "ACGTACGACGTACGACGTACG";
	
	make_seq(ss.begin(),ss.end(),tree);
	
	for(int i=0;i<ss.length();++i) {
		ul_find[i] = i;
		d_find[i] = i;
		w_find[i] = dawg::evo_node_weight<>(i,i);
	}
			
	FT tree2(tree);
	FT tree3;
	tree3 = tree2;
	
	printnode<dawg::residue, dawg::evo_node_weight<> >(&(*tree.root()));
	cout << endl;
	printnode<dawg::residue, dawg::evo_node_weight<> >(&(*tree2.root()));
	cout << endl;
	printnode<dawg::residue, dawg::evo_node_weight<> >(&(*tree3.root()));
	cout << endl;
		
	for(int i=0;i<ss.length();++i) {
		cout << make_seq.decode(tree[ul_find[i]].base());
	}
	cout << endl;

	for(int i=0;i<ss.length();++i) {
		cout << make_seq.decode(tree2[ul_find[i]].base());
	}
	cout << endl;

	for(int i=0;i<ss.length();++i) {
		cout << make_seq.decode(tree3[ul_find[i]].base());
	}
	cout << endl;

	
	return 0;
}

/*
	srand(145);
	for(int i=0;i<10;++i) {
		store.push_back(dawg::residue(rand()&3,1.0,0,1.0));
		cout << store.back().base();
	}
	cout << endl;
	tree.insert(tree.begin(),store.begin(),store.end());
	
	for(tit=tree.begin();tit!=tree.end();++tit)
		cout << make_seq.decode(tit->val.base());
	cout << endl;
*/

