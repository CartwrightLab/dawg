#pragma once

#include <vector>
#include <string>
#include <algorithm>

unsigned long hash_adler32(const char* cs);

class Node {
public:
	Node(const char *cs, double d, Node* p)
	{
		dLen = d;
		pSub = p;
		pSib = NULL;
		if(cs)
			strcpy(csLabel, cs);
		else
		{
			std::vector<std::string> v;
			while(p != NULL)
			{
				v.push_back(p->csLabel);
				p = p->pSib;
			}
			std::sort(v.begin(), v.end());
			std::string ss;
			for(std::vector<std::string>::iterator it = v.begin();
				it != v.end(); ++it)
			{
				ss += *it;
				ss += ",";
			}
			csLabel[0] = '!';
			_ltoa(hash_adler32(ss.c_str()), csLabel+1, 16);
		}
	}
	~Node()
	{
		if(pSub)
			delete pSub;
		if(pSib)
			delete pSib;
	}
	
	char csLabel[1024];
	double dLen;
	Node* pSub;
	Node* pSib;
	static std::vector<Node*> s_stack;
};