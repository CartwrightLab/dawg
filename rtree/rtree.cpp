// rtree.cpp : Defines the entry point for the console application.
//

#include <tchar.h>
#include "rtree.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <iomanip>

using namespace std;

vector<Node*> Node::s_stack;

unsigned long hash_adler32(const char* cs);

bool Parse(const char* cs);

class RTree
{
public:
	struct Section;
	typedef	vector<Section> RNode;
	typedef map<string, RNode> Map;
	struct Section
	{
		Section() { }
		Section(double d, Map::iterator it) : dBranchLen(d), itAncestor(it) { }
		double dBranchLen;
		Map::iterator itAncestor;
	};
	
	RTree() : m_nSec(0) {}

	Map m_map;
	int m_nSec;
	void Process(Node* pNode)
	{
		ProcessNode(pNode);
		m_nSec++;
	}

	void ProcessNode(Node* pNode)
	{
		static vector<Map::iterator> vStack;

		if(pNode->pSib)
			ProcessNode(pNode->pSib);

		Section sec;
		sec.dBranchLen = pNode->dLen;
		if(vStack.empty())
			sec.itAncestor = m_map.end();
		else
			sec.itAncestor = vStack.back();
		RNode& rn = m_map[pNode->csLabel];
		for(RNode::const_iterator it = rn.begin(); it != rn.end(); ++it)
		{
			if(sec.itAncestor == it->itAncestor)
				sec.dBranchLen = it->dBranchLen;
		}
		rn.resize(m_nSec+1, Section(0.0, m_map.end()));
		rn[m_nSec] = sec;
				
		if(pNode->pSub)
		{
			vStack.push_back(m_map.find(pNode->csLabel));
			ProcessNode(pNode->pSub);
			vStack.pop_back();
		}
	}
};

//class RTree
//{
//public:
//	vector<RNode> m_vNodes;
//	vector<string> m_vLabels;
//
//	void Process(Node* pNode)
//	{
//		static vector<int> vStack;
//
//		if(pNode->pSib)
//			Process(pNode->pSib);
//
//		Section sec;
//		sec.dBranchLen = pNode->dLen;
//		if(vStack.empty())
//			sec.nAncestor = -1;
//		else
//			sec.nAncestor = vStack.back();
//		vector<string>::iterator pos = find(m_vLabels.begin(), m_vLabels.end(), pNode->csLabel);
//		if(pos == m_vLabels.end())
//		{
//			m_vLabels.push_back(pNode->csLabel);
//			m_vNodes.push_back(RNode());
//			pos = m_vLabels.end()-1;
//		}
//		int id = pos-m_vLabels.begin();
//		m_vNodes[id].push_back(sec);
//		
//		if(pNode->pSub)
//		{
//			vStack.push_back(id);
//			Process(pNode->pSub);
//			vStack.pop_back();
//		}
//	}
//};

int _tmain(int argc, _TCHAR* argv[])
{
	Parse("sample.txt");
	RTree tree;
	for(vector<Node*>::iterator it = Node::s_stack.begin();
		it != Node::s_stack.end(); ++it)
	{
		tree.Process(*it);
	}
	for(RTree::Map::iterator it = tree.m_map.begin(); it != tree.m_map.end(); ++it)
	{
		cout << it->first << endl;
		for(int j=0; j< it->second.size(); ++j)
		{
			cout << "    " << j << ": " << it->second[j].dBranchLen << " ";
			if(it->second[j].itAncestor == tree.m_map.end())
				cout << "!NONE" << endl;
			else
				cout << it->second[j].itAncestor->first << endl;

		}
	}
	//for(int i=0; i<tree.m_vNodes.size(); ++i)
	//{
	//	cout << tree.m_vLabels[i] << " " << setbase(16) << hash_adler32(tree.m_vLabels[i].c_str()) <<  endl;
	//	for(int j=0; j< tree.m_vNodes[i].size(); ++j)
	//	{
	//		cout << "    " << j << ": ";
	//		if(tree.m_vNodes[i][j].nAncestor == -1)
	//			cout << "NULL" << endl;
	//		else
	//			cout << tree.m_vLabels[tree.m_vNodes[i][j].nAncestor] << endl;
	//	}
	//}

	
	for(vector<Node*>::iterator it = Node::s_stack.begin();
		it != Node::s_stack.end(); ++it)
		delete *it;
	return 0;
}

unsigned long hash_adler32(const char* cs)
{
	const unsigned long BASE = 65521;
	unsigned long s1 = 1;
	unsigned long s2 = 0;
	size_t len = strlen(cs);

	if (len % 8 != 0)
	{
		do
		{
			s1 += *cs++;
			s2 += s1;
			len--;
		} while (len % 8 != 0);

		if (s1 >= BASE)
			s1 -= BASE;
		s2 %= BASE;
	}
	while (len > 0)
	{
		s1 += cs[0]; s2 += s1;
		s1 += cs[1]; s2 += s1;
		s1 += cs[2]; s2 += s1;
		s1 += cs[3]; s2 += s1;
		s1 += cs[4]; s2 += s1;
		s1 += cs[5]; s2 += s1;
		s1 += cs[6]; s2 += s1;
		s1 += cs[7]; s2 += s1;

		len -= 8;
		cs += 8;

		if (s1 >= BASE)
			s1 -= BASE;
		if (len % 0x8000 == 0)
			s2 %= BASE;
	}
	return (s1 << 16) | (s2 & 0xFFFF);
}

