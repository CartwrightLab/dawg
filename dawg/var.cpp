#include "dawg.h"
#include <map>
#include <float.h>
#include <stdio.h>

using namespace std;

DawgVar::~DawgVar()
{
	Unset();
}

class MapSsToVar : public std::map<std::string, DawgVar*>
{
public:	
	virtual ~MapSsToVar()
	{
		for(iterator it = begin(); it != end(); ++it)
			delete it->second;
		clear();
	}
};

void DawgVar::ClearMap()
{
	for(MapSsToVar::iterator it = GetMap().begin(); it != GetMap().end(); ++it)
			delete it->second;
	GetMap().clear();
}

MapSsToVar& DawgVar::GetMap()
{
	static MapSsToVar s_map;
	return s_map;
}

DawgVar* DawgVar::GetVar(const std::string &ssKey)
{
	MapSsToVar::iterator it = GetMap().find(ssKey);
	return (it != GetMap().end()) ? it->second : NULL;
}

void DawgVar::SetVar(const std::string &ssKey, DawgVar* pVar)
{
	GetMap()[ssKey] = pVar;
}

//bool DawgVar::Get(const std::string& ssKey, double &rdVar )
//{
//	DawgVar *pVar = GetVar(ssKey);
//	return ( pVar != NULL && pVar->Get(rdVar));
//}
//
//bool DawgVar::Get(const std::string& ssKey, int &rnVar )
//{
//	DawgVar *pVar = GetVar(ssKey);
//	return ( pVar != NULL && pVar->Get(rnVar));
//}
//
//bool DawgVar::Get(const std::string& ssKey, std::string &rssVar )
//{
//	DawgVar *pVar = GetVar(ssKey);
//	return ( pVar != NULL && pVar->Get(rssVar));
//}
//
//bool DawgVar::Get(const std::string& ssKey, bool &rbVar )
//{
//	DawgVar *pVar = GetVar(ssKey);
//	return ( pVar != NULL && pVar->Get(rbVar));
//}
//
//bool DawgVar::Get(const std::string &ssKey, Vec *&rpVec)
//{
//	DawgVar *pVar = GetVar(ssKey);
//	return ( pVar != NULL && pVar->Get(rpVec));
//}

bool DawgVar::Get( double &rdVar ) const
{
	if(!IsType(tyNumber))
		return false;
	rdVar = GetNumber();
	return true;
}

bool DawgVar::Get( int &rnVar ) const
{
	if(!IsType(tyNumber))
		return false;
	rnVar = (int)GetNumber();
	return true;
}

bool DawgVar::Get(bool &rbVar ) const
{
	if(!IsType(tyBool))
		return false;
	rbVar = GetBool();
	return true;
}

bool DawgVar::Get(std::string &rssVar ) const
{
	if(!IsType(tyString))
		return false;
	rssVar = GetString();
	return true;
}

bool DawgVar::Get(const Vec *&rVec) const
{
	if(!IsType(tyVector))
		return false;
	rVec = GetVector();
	return true;
}

void DawgVar::Unset()
{
	switch(GetType())
	{
		case tyString:
			delete m_pssData; break;
		case tyVector:
			delete m_pvData;  break;
		case tyTree:
			delete m_ptrData; break;
	}
	m_tyType = tyNone;
}

void DawgVar::Set(double dVar )
{
	Unset();
	m_tyType = tyNumber;
	m_dData = dVar;
}

void DawgVar::Set(int nVar )
{
	Unset();
	m_tyType = tyNumber;
	m_dData = (double)nVar;
}

void DawgVar::Set(bool bVar )
{
	Unset();
	m_tyType = tyBool;
	m_bData = bVar;
}

void DawgVar::Set(const std::string &rssVar )
{
	Unset();
	m_tyType = tyString;
	m_pssData = new string(rssVar);
}

void DawgVar::Set(const Vec *pVec)
{
	Unset();
	m_tyType = tyVector;
	m_pvData = pVec;
}

bool DawgVar::Get(Tree *&rTree) const
{
	if(!IsType(tyTree))
		return false;
	rTree = GetTree();
	return true;
}
void DawgVar::Set(Tree *pTree)
{
	Unset();
	m_tyType = tyTree;
	m_ptrData = pTree;
}
//bool DawgVar::Get(const string & ssKey, Tree *&rTree)
//{
//	DawgVar *pVar = GetVar(ssKey);
//	return ( pVar != NULL && pVar->Get(rTree));
//}

unsigned int DawgVar::Size()
{
	switch(m_tyType)
	{
		case tyVector:  return m_pvData->size();
		case tyNone:	return 0;
		default:		return 1;
	};
}


/********************************************************
class Node
********************************************************/
Node::~Node()
{
	if(m_pSib != NULL)
		delete m_pSib;
	if(m_pChild != NULL)
		delete m_pChild;
}

string Node::ToString() const
{
	string ssTemp = "";
	if(IsClade())
	{
		ssTemp += '(';
		ssTemp += m_pChild->ToString();
		Node* p = m_pChild->m_pSib;
		while(p)
		{
			ssTemp += ',';
			ssTemp += p->ToString();
			p = p->m_pSib;
		}
		ssTemp += ')';
	}
	ssTemp += Label();
	double d = BranchLength();
	if(d != 0.0)
	{
		ssTemp += ':';
		char csBuffer[32];
		sprintf(csBuffer, "%0.10f", BranchLength());
		char *p = &csBuffer[0];
		while (*p) {p++; }
		while(*(--p) == '0' || *p == '.')
			*p = '\0';
		ssTemp += csBuffer;
	}
	return ssTemp;
}
