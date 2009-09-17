#include <sys/time.h>

#include <iostream>
#include <string>


#include "nseq.h"

using namespace std;

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
extern FT::iterator tit;
FT::iterator tit;

int main(int argc, char* argv[]) {
	FT tree;
	char in;
	dawg::residue_factory make_seq;
	
	std::string ss = "acgtacgtacgtACGTACGTaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	
	make_seq(ss.begin(),ss.end(),tree);
	
	for(int i=0;i<ss.length();++i) {
		ul_find[i] = i;
		d_find[i] = i;
		w_find[i] = dawg::evo_node_weight<>(i,i);
	}
		
	for(int i=0;i<ss.length();++i) {;
		cout << make_seq.decode(tree[ul_find[i]].base());
//		cout << tit->val.base() << " ";
//		tit = tree.find(d_find[i]);
//		cout << tit->val.base() << " ";
//		tit = tree.find(w_find[i]);
//		cout << tit->val.base() << endl;
	}
	cout << endl;
		
	return 0;
}

