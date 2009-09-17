#include <iostream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>


#include "nseq.h"

using namespace std;
using namespace boost::posix_time;

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
extern ptime tstart,tend;
ptime tstart,tend;

vector<dawg::residue> store;

int main(int argc, char* argv[]) {
	FT tree;
	char in;
	dawg::residue_factory make_seq;
	
	std::string ss = "ACGTACGTACGTACGTACGTACGTACGTACGT";
	
	make_seq(ss.begin(),ss.end(),tree);
	
	for(int i=0;i<ss.length();++i) {
		ul_find[i] = i;
		d_find[i] = i;
		w_find[i] = dawg::evo_node_weight<>(i,i);
	}
		
	for(int i=0;i<ss.length();++i) {
		cout << make_seq.decode(tree[ul_find[i]].base());
//		cout << tit->val.base() << " ";
//		tit = tree.find(d_find[i]);
//		cout << tit->val.base() << " ";
//		tit = tree.find(w_find[i]);
//		cout << tit->val.base() << endl;
	}
	cout << endl;
	
	for(int i=0;i<1000;++i) {
		store.push_back(dawg::residue(3,1.0,0,1.0));
	}
	store.push_back(dawg::residue(2,1.0,0,1.0));
	
	tree.insert(tree.begin(),store.begin(),store.end());
	tree.insert2(tree.end(),store.begin(),store.end());
	
	for(tit=tree.begin();tit!=tree.end();++tit)
		cout << make_seq.decode(tit->val.base());
	cout << endl;
		
	tstart = microsec_clock::local_time();
	tend = microsec_clock::local_time();
	cout << "Time 0: " << (tend-tstart).total_microseconds()/1e6 << endl;
	
	for(int ii=0;ii<5;++ii) {
	tstart = microsec_clock::local_time();
	for(int i=0; i<100;++i) {
		FT t;
		make_seq(ss.begin(),ss.end(),t);
		tit = t.find(10ul);
		for(int j=0; j<500;++j) {
			t.insert2(tit, store.begin(), store.end());
		}
	}
	tend = microsec_clock::local_time();
	cout << "Time 1: " << (tend-tstart).total_microseconds()/1e6 << endl;

	tstart = microsec_clock::local_time();
	for(int i=0; i<100;++i) {
		FT t;
		make_seq(ss.begin(),ss.end(),t);
		tit = t.find(10ul);
		for(int j=0; j<500;++j) {
			t.insert(tit, store.begin(), store.end());
		}
	}
	tend = microsec_clock::local_time();
	cout << "Time 2: " << (tend-tstart).total_microseconds()/1e6 << endl;
	}
	return 0;
}

