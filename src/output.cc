// output.cc - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"

using namespace std;

void PrintVar(const DawgVar* var);

void   PrintSequencesFasta(ostream &os, const Tree::Alignment& aln);
void   PrintSequencesNexus(ostream &os, const Tree::Alignment& aln);
void   PrintSequencesPhylip(ostream &os, const Tree::Alignment& aln);
void   PrintSequencesClustal(ostream &os, const Tree::Alignment& aln);

FileFormat g_fileFormat = FASTA;
const char *g_csBlock = NULL;
int			g_nDataSet = 0;
int			g_nDataSetNum = 1;

bool SetFormat(FileFormat fmt, int nNum, const char* csBlock)
{
	g_fileFormat = fmt;
	g_csBlock = csBlock;
	g_nDataSet = 0;
	g_nDataSetNum = nNum;
	return true;
}

void DawgIniOutput(ostream& os)
{
	switch(g_fileFormat)
	{
	case NEXUS:
		os << "#NEXUS" << endl << "[Created by DAWG Version "
			<< VERSION << endl << endl;
		os << /*EvoDescription() <<*/ ']' << endl;
		break;
	case CLUSTAL:
		os << "CLUSTAL multiple sequence alignment (Created by DAWG Version "
			<< VERSION << ")" << endl << endl << endl;
		break;
	case FASTA:
	case PHYLIP:
	default:
		break;
	}

}

bool SaveAlignment(ostream &rFile, const Tree::Alignment& aln)
{
	if(g_nDataSetNum > 1)
		rFile << "[DataSet " << ++g_nDataSet << ']' << endl;
	switch(g_fileFormat)
	{
		case NEXUS:
			PrintSequencesNexus(rFile, aln);
			break;
		case PHYLIP:
			PrintSequencesPhylip(rFile, aln);
			break;
		case CLUSTAL:
			PrintSequencesClustal(rFile, aln);
			break;
		case FASTA:
		default:
			PrintSequencesFasta(rFile, aln);
	};
	return true;
}

void PrintSequencesFasta(ostream& os, const Tree::Alignment& aln)
{
	for(Tree::Alignment::const_iterator cit = aln.begin(); cit != aln.end(); ++cit)
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

void PrintSequencesNexus(ostream &os, const Tree::Alignment& aln)
{
	Tree::Alignment::const_iterator cit = aln.begin();
	os << "BEGIN DATA;" << endl;
	os << "\tDIMENSIONS NTAX=" << aln.size() << " NCHAR=" << cit->second.length() << ';' << endl;
	os << "\tFORMAT MISSING=? GAP=- DATATYPE=DNA;" << endl;
	os << "\tMATRIX" << endl;
	for(; cit != aln.end(); ++cit)
		os << setw(15) << setiosflags(ios::left) << cit->first << ' ' << cit->second << endl;
	os << ';' << endl << "END;" << endl << endl;
	if(g_csBlock)
		os << g_csBlock << endl << endl;
}
void PrintSequencesPhylip(ostream &os, const Tree::Alignment& aln)
{
	Tree::Alignment::const_iterator cit = aln.begin();
	os << aln.size() << ' ' << cit->second.length() << endl;
	os << setfill(' ');
	for(; cit != aln.end(); ++cit)
		os << setw(10) << setiosflags(ios::left) << cit->first.substr(0,10) << cit->second << endl;
	os << endl;
}

void   PrintSequencesClustal(ostream &os, const Tree::Alignment& aln)
{
	unsigned long uLen = aln.begin()->second.length();
	unsigned long l;
	for(unsigned long u = 0; u < uLen; u+=l)
	{
		l = min(60ul, uLen);
		for(Tree::Alignment::const_iterator cit = aln.begin(); cit != aln.end(); ++cit)
			os << setw(15) << setiosflags(ios::left) << cit->first << " " << cit->second.substr(u, l) << endl;
		os << endl << endl;
	}
}
