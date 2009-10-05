#include <iostream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>


#include "sequence.h"

using namespace std;

typedef dawg::finger_tree<dawg::residue, dawg::evo_node_weigher<> > FT;

void printnode(FT::node::pointer p) {
	if(p == NULL)
		return;
	cout << "(";
	if(p->left != p)
		printnode(p->left);
	cout << p->val.base() << ((p->color) ? "r" : "b") << (p->weight.length);
	if(p->right != p)
		printnode(p->right);
	cout << ")";
}



int main(int argc, char* argv[]) {
	FT tree;
	char in;
	dawg::residue_factory make_seq;

	std::string ss = "ACGTACGTACGT";

	for(std::string::iterator it = ss.begin(); it != ss.end(); ++it) {
		tree.insert(tree.end(), make_seq(*it));
	}

	cout << tree.root()->weight.length << endl;
	printnode(&*tree.root());
	cout << endl;

	FT::iterator it = tree.find(FT::weight_type::size_type(4));
	unsigned int u = 8;
	while(u--) {
		it = tree.find(FT::weight_type::size_type(u));
		//it->val.mark_deleted(true);
		cout << u << " " << it->val.base() <<  " " << tree.root()->weight.length << endl;
		printnode(&*tree.root());
		cout << endl;
		//it = tree.search_and_update(it, FT::size_type(0));
	}
	//it->update_weight();
	cout << it->val.base() <<  " " << tree.root()->weight.length << endl;
	printnode(&*tree.root());
	cout << endl;

	return 0;
}

