// var.cc - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"

using namespace std;

DawgVar::~DawgVar()
{
	Unset();
}

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

bool DawgVar::Get(NewickNode *&rTree) const
{
	if(!IsType(tyTree))
		return false;
	rTree = GetTree();
	return true;
}
void DawgVar::Set(NewickNode *pTree)
{
	Unset();
	m_tyType = tyTree;
	m_ptrData = pTree;
}

DawgVar::Vec::size_type DawgVar::Size()
{
	switch(m_tyType)
	{
		case tyVector:  return m_pvData->size();
		case tyNone:	return 0;
		default:		return 1;
	};
}
