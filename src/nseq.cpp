#include <iostream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>


#include "sequence.h"

using namespace std;

typedef dawg::finger_tree<dawg::residue, dawg::evo_node_weight<> > FT;

void printnode(FT::node::pointer p) {
	if(p == NULL)
		return;
	cout << "(";
	printnode(p->left);
	cout << p->val.base() << ((p->color) ? "r" : "b") << (p->weight.length);
	printnode(p->right);
	cout << ")";
}



int main(int argc, char* argv[]) {
	FT tree;
	char in;
	dawg::residue_factory make_seq;

	std::string ss = "ACGTACGTACGT";

	make_seq(ss.begin(),ss.end(),tree);

	cout << tree.root()->weight.length << endl;
	printnode(&*tree.root());
	cout << endl;

	FT::iterator it = tree.find(4ul);
	unsigned int u = 8;
	while(u--) {
		it->val.mark_deleted(true);
		cout << it->val.base() <<  " " << tree.root()->weight.length << endl;
		printnode(&*tree.root());
		cout << endl;
		it = tree.search_and_update(it, 0ul);
	}
	//it->update_weight();
	cout << it->val.base() <<  " " << tree.root()->weight.length << endl;
	printnode(&*tree.root());
	cout << endl;

	return 0;
}

