#include "dawg.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <assert.h>

using namespace std;

void PrintVar(const DawgVar* var);
void PrintSeq(const Seq& rSeq);
void PrintTreeSeq(Tree* pTree);
void PrintTree(const Node* pNode);

string SeqToString(const Seq& rSeq);

typedef map<string, string> MapLabelToAln;
typedef map<const Node*, string> MapNodeToSs;

void   FillSequences(MapNodeToSs& rMap, const Tree* pTree);
void   MapSequences(MapNodeToSs& rMap, const Node* pNode);
void   MapAlignment(MapLabelToAln& rMap, const Tree* pTree);
void   Align(MapNodeToSs &rMap);
void   PrintSequencesFasta(ostream &os, const MapLabelToAln& rMap);
void   PrintSequencesNexus(ostream &os, const MapLabelToAln& rMap);
void   PrintSequencesPhylip(ostream &os, const MapLabelToAln& rMap);

FileFormat g_fileFormat = FASTA;
const char *g_csBlock = NULL;
int			g_nDataSet = 0;
int			g_nDataSetNum = 1;
bool		g_bGapSingle = false;

bool SetFormat(FileFormat fmt, int nNum, const char* csBlock, bool bGapSingle)
{
	g_fileFormat = fmt;
	g_csBlock = csBlock;
	g_nDataSet = 0;
	g_bGapSingle = bGapSingle;
	g_nDataSetNum = nNum;
	return true;
}

bool DawgOpen(const char* csFile, ofstream& rFile)
{
	rFile.open(csFile);
	if(!rFile.is_open())
		return false;
	return true;
}

void DawgIniOutput(ostream& os)
{
	if(g_fileFormat == NEXUS)
	{
		os << "#NEXUS" << endl << "[Created by DAWG Version " << DawgVersion() << endl << endl;
		os << EvoDescription() << ']' << endl;
	}
}

bool SaveSequences(ostream &rFile, const Tree* arTrees[], unsigned int uSize)
{
	++g_nDataSet;
	MapLabelToAln m;
	for(unsigned int u=0;u<uSize;++u)
		MapAlignment(m, arTrees[u]);
	switch(g_fileFormat)
	{
		case FASTA:
			PrintSequencesFasta(rFile, m);
			break;
		case NEXUS:
			PrintSequencesNexus(rFile, m);
			break;
		case PHYLIP:
			PrintSequencesPhylip(rFile, m);
			break;
		default:
			PrintSequencesFasta(rFile, m);
	};
	return true;
}

string SeqToString(const Seq& rSeq)
{
	string str;
	str.reserve(rSeq.size());
	for(Seq::const_iterator nit = rSeq.begin(); nit != rSeq.end(); ++nit)
		str += NucToChar(nit->m_nuc);
	return str;
}

// Appends aligned sequeces from pTree to the map
void MapAlignment(MapLabelToAln& rMap, const Tree* pTree)
{
	//Construct map of template strings
	MapNodeToSs ssMap;
	MapSequences(ssMap, pTree);

	//Align Sequences
    Align(ssMap);
	FillSequences(ssMap, pTree);

	//Append Sequences with label in rMap
	for(MapNodeToSs::iterator it = ssMap.begin(); it != ssMap.end(); ++it)
		if(!it->first->Label().empty())
			rMap[it->first->Label()] += it->second;
}

void MapSequences(MapNodeToSs& rMap, const Node* pNode)
{
	rMap[pNode] = pNode->Gaps();
	if(pNode->Child())
		MapSequences(rMap, pNode->Child());
	if(pNode->Sibling())
		MapSequences(rMap, pNode->Sibling());
}

// Insertion: ^
// Deletion:  -
// Deleted Insertion: =
// Root Nuc: %

void Align(MapNodeToSs &rMap)
{
	unsigned int u=0;
	bool bEnd=true;
	while(bEnd)
	{
		bool bIns=false, bBeyond=false;
		bEnd = false;
		for(MapNodeToSs::iterator it=rMap.begin();it!=rMap.end();++it)
		{
			if(u >= it->second.length())
			{
				//bBeyond=true;
				continue;
			}
			bEnd = true;
			if(it->second[u] == '^' || it->second[u] == '=')
			{
				bIns = true;
				break;
			}
		}
		if(bIns || (bBeyond && bEnd) )
		{
			for(MapNodeToSs::iterator it=rMap.begin();it!=rMap.end();++it)
			{
				if( u >= it->second.length())
					it->second.resize(u+1, '-');
				else if(it->second[u] != '^' && it->second[u] != '=')
					it->second.insert(u, 1, '-');
			}
		}
		u++;
	}
}

// Convert sequences to human format and store in map under node.
void FillSequences(MapNodeToSs& rMap, const Tree* pTree)
{
	string ssTemp = SeqToString(pTree->Sequence());
	string& rss = rMap[pTree];
	unsigned int v=0;
	char bPrevChar = '\0';
	for(string::iterator it=rss.begin(); it != rss.end(); ++ it)
	{
		if(*it == '%' || *it == '^')
			*it = ssTemp[v++];
		else if(*it == '=' || *it == '-')
			*it = (g_bGapSingle && (bPrevChar == '-' || bPrevChar == '?')) ? '?' : '-';
		bPrevChar = *it;
	}
	if(pTree->Child())
		FillSequences(rMap, pTree->Child());
	if(pTree->Sibling())
		FillSequences(rMap, pTree->Sibling());
}

void PrintSequencesFasta(ostream& os, const MapLabelToAln& rMap)
{
	if(g_nDataSetNum > 1)
		os << "[DataSet " << g_nDataSet << ']' << endl;
	for(MapLabelToAln::const_iterator cit = rMap.begin(); cit != rMap.end(); ++cit)
	{
		// Print Label
		os << '>' << cit->first << endl;
		// Print Sequence
		unsigned int u=0;
		while(u < cit->second.length())
		{
			os << cit->second[u];
			if((++u % 60) == 0)
				os << endl;
		}
		if(u%60)
			os << endl;
		os << endl;
	}
}

void PrintSequencesNexus(ostream &os, const MapLabelToAln& rMap)
{
	MapLabelToAln::const_iterator cit = rMap.begin();	
	os << "BEGIN DATA; [DataSet " << g_nDataSet << ']' << endl;
	os << "\tDIMENSIONS NTAX=" << rMap.size() << " NCHAR=" << cit->second.length() << ';' << endl;
	os << "\tFORMAT MISSING=? GAP=- DATATYPE=DNA;" << endl;
	os << "\tMATRIX" << endl;
	for(; cit != rMap.end(); ++cit)
		os << setw(15) << setiosflags(ios::left) << cit->first << ' ' << cit->second << endl;
	os << ';' << endl << "END;" << endl << endl;
	if(g_csBlock)
		os << g_csBlock << endl << endl;
}
void PrintSequencesPhylip(ostream &os, const MapLabelToAln& rMap)
{
	if(g_nDataSetNum > 1)
		os << "[DataSet " << g_nDataSet << ']' << endl;
	MapLabelToAln::const_iterator cit = rMap.begin();
	os << rMap.size() << ' ' << cit->second.length() << endl;
	os << setfill(' ');
	for(; cit != rMap.end(); ++cit)
		os << setw(10) << setiosflags(ios::left) << cit->first.substr(0,10) << cit->second << endl;
	os << endl;
}

void PrintTree(const Node* pNode)
{
	cout << pNode->ToString();
}

void PrintVar(const DawgVar* var)
{
	switch(var->GetType())
	{
	case DawgVar::tyNone:
		cout << "(null)";
		break;
	case DawgVar::tyString:
		cout << var->GetString();
		break;
	case DawgVar::tyNumber:
		cout << var->GetNumber();
		break;
	case DawgVar::tyBool:
		cout << var->GetBool();
		break;
	case DawgVar::tyVector:
		{
			cout << '{';
			const DawgVar::Vec *v = var->GetVector();
			if(!v->empty())
			{
				DawgVar::Vec::const_iterator it = v->begin();
				PrintVar(*it);
				for( ++it;it != v->end(); ++it)
				{
					cout << ",";
					PrintVar(*it);
				}
			}
			cout << '}';
			break;
		}
	case DawgVar::tyTree:
		{
			PrintTree(var->GetTree());
			cout << ";";
			break;
		}
	}
}

