#include "dawg.h"
#include "var.h"

using namespace std;

void PrintVar(const DawgVar* var);

void   PrintSequencesFasta(ostream &os, const Tree::Alignment& aln);
void   PrintSequencesNexus(ostream &os, const Tree::Alignment& aln);
void   PrintSequencesPhylip(ostream &os, const Tree::Alignment& aln);

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
		os << "#NEXUS" << endl << "[Created by DAWG Version " << VERSION << endl << endl;
		os << EvoDescription() << ']' << endl;
	}
}

bool SaveAlignment(ostream &rFile, const Tree::Alignment& aln)
{
	++g_nDataSet;
	switch(g_fileFormat)
	{
		case FASTA:
			PrintSequencesFasta(rFile, aln);
			break;
		case NEXUS:
			PrintSequencesNexus(rFile, aln);
			break;
		case PHYLIP:
			PrintSequencesPhylip(rFile, aln);
			break;
		default:
			PrintSequencesFasta(rFile, aln);
	};
	return true;
}

void PrintSequencesFasta(ostream& os, const Tree::Alignment& aln)
{
	if(g_nDataSetNum > 1)
		os << "[DataSet " << g_nDataSet << ']' << endl;
	for(Tree::Alignment::const_iterator cit = aln.begin(); cit != aln.end(); ++cit)
	{
		// Skip any sequence that begin with one of the two special characters
		if(cit->first[0] == '(' || cit->first[0] == '_')
			continue;
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

void PrintSequencesNexus(ostream &os, const Tree::Alignment& aln)
{
	Tree::Alignment::const_iterator cit = aln.begin();	
	os << "BEGIN DATA; [DataSet " << g_nDataSet << ']' << endl;
	os << "\tDIMENSIONS NTAX=" << aln.size() << " NCHAR=" << cit->second.length() << ';' << endl;
	os << "\tFORMAT MISSING=? GAP=- DATATYPE=DNA;" << endl;
	os << "\tMATRIX" << endl;
	for(; cit != aln.end(); ++cit)
	{
		if(cit->first[0] == '(' || cit->first[0] == '_')
			continue;
		os << setw(15) << setiosflags(ios::left) << cit->first << ' ' << cit->second << endl;
	}
	os << ';' << endl << "END;" << endl << endl;
	if(g_csBlock)
		os << g_csBlock << endl << endl;
}
void PrintSequencesPhylip(ostream &os, const Tree::Alignment& aln)
{
	if(g_nDataSetNum > 1)
		os << "[DataSet " << g_nDataSet << ']' << endl;
	Tree::Alignment::const_iterator cit = aln.begin();
	os << aln.size() << ' ' << cit->second.length() << endl;
	os << setfill(' ');
	for(; cit != aln.end(); ++cit)
	{
		if(cit->first[0] == '(' || cit->first[0] == '_')
			continue;
		os << setw(10) << setiosflags(ios::left) << cit->first.substr(0,10) << cit->second << endl;
	}
	os << endl;
}

//void PrintTree(const Node* pNode)
//{
//	cout << pNode->ToString();
//}
//
//void PrintVar(const DawgVar* var)
//{
//	switch(var->GetType())
//	{
//	case DawgVar::tyNone:
//		cout << "(null)";
//		break;
//	case DawgVar::tyString:
//		cout << var->GetString();
//		break;
//	case DawgVar::tyNumber:
//		cout << var->GetNumber();
//		break;
//	case DawgVar::tyBool:
//		cout << var->GetBool();
//		break;
//	case DawgVar::tyVector:
//		{
//			cout << '{';
//			const DawgVar::Vec *v = var->GetVector();
//			if(!v->empty())
//			{
//				DawgVar::Vec::const_iterator it = v->begin();
//				PrintVar(*it);
//				for( ++it;it != v->end(); ++it)
//				{
//					cout << ",";
//					PrintVar(*it);
//				}
//			}
//			cout << '}';
//			break;
//		}
//	case DawgVar::tyTree:
//		{
//			PrintTree(var->GetTree());
//			cout << ";";
//			break;
//		}
//	}
//}

