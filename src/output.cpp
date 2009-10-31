// output.cc - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"

using namespace std;

void FilterAlignment(const Tree::Alignment& alnSrc, Tree::Alignment& alnDest, unsigned int uFlags);
void FilterSequence(const string& ssSrc, string& ssDest, unsigned int uFlags);
void PrintSequencesFasta(  ostream &os, const Tree::Alignment& aln, unsigned int uFlags);
void PrintSequencesNexus(  ostream &os, const Tree::Alignment& aln, unsigned int uFlags);
void PrintSequencesPhylip( ostream &os, const Tree::Alignment& aln, unsigned int uFlags);
void PrintSequencesClustal(ostream &os, const Tree::Alignment& aln, unsigned int uFlags);

unsigned int g_fileFormat = FormatFasta;
const char *g_csHead = NULL;
const char *g_csBefore = NULL;
const char *g_csAfter = NULL;
const char *g_csTail = NULL;
int			g_nDataSet = 0;
int			g_nDataSetNum = 1;
bool		g_bSubst = true;

// Reset Format
bool SetFormat(unsigned int fmt, int nNum,
			   const char* csHead, const char* csBefore,
			   const char* csAfter, const char* csTail,
			   bool bSubst)
{
	g_fileFormat = fmt;
	g_csHead = csHead;
	g_csBefore = csBefore;
	g_csAfter = csAfter;
	g_csTail = csTail;
	g_nDataSet = 0;
	g_nDataSetNum = nNum;
	g_bSubst = bSubst;
	return true;
}

string Subst(const char* cs)
{
	if(!g_bSubst)
		return cs;
	string ss = "";
	char csBuffer[65] = "";
	const char *cs1 = cs;
	const char *cs2 = strchr(cs1, '%');
	while(cs2 != NULL)
	{
		ss.append(cs1, cs2-cs1);
		switch(cs2[1])
		{
		case 'r':
			snprintf(csBuffer, 64, "%d", g_nDataSet);
			ss.append(csBuffer);
			cs1 = cs2+2;
			break;
		case 'R':
			snprintf(csBuffer, 64, "%d", g_nDataSetNum);
			ss.append(csBuffer);
			cs1 = cs2+2;
			break;
		case '%':
			ss.append("%");
			cs1 = cs2+2;
			break;
		default:
			ss.append("%");
			cs1 = cs2+1;
			break;
		}

		cs2 = strchr(cs1, '%');
	}
	ss.append(cs1);
	return ss;
}

// Initialize Output for formats that need it
void DawgIniOutput(ostream& os)
{
	switch(g_fileFormat)
	{
	case FormatNexus:
		os << "#NEXUS" << endl << "[Created by " << PACKAGE_STRING << ']' << endl;
		break;
	case FormatClustal:
		os << "CLUSTAL multiple sequence alignment (Created by " << PACKAGE_STRING << ")" << endl << endl << endl;
		break;
	}
	if(g_csHead != NULL)
		os << Subst(g_csHead) << endl;
}

// Finalize Output
void DawgFinOutput(ostream& os)
{
	if(g_csTail != NULL)
		os << Subst(g_csTail) << endl;
}


// Save alignment to output stream
bool SaveAlignment(ostream &rFile, const Tree::Alignment& aln, unsigned int uFlags)
{
	g_nDataSet++;
	if(g_csBefore != NULL)
		rFile << Subst(g_csBefore) << endl;
	// Filter sequences to a local alignment, applying format modifiers
	Tree::Alignment alnLocal;
	FilterAlignment(aln, alnLocal, uFlags);
	for(Tree::Alignment::const_iterator cit = alnLocal.begin(); cit != alnLocal.end(); ++cit)
		FilterSequence(cit->second, alnLocal[cit->first], uFlags);
	
	// Save Alignment in requested file format
	// Could be an array if flexibility is needed
	switch(g_fileFormat)
	{
		case FormatNexus:
			PrintSequencesNexus(rFile, alnLocal, uFlags);
			break;
		case FormatPhylip:
			PrintSequencesPhylip(rFile, alnLocal, uFlags);
			break;
		case FormatClustal:
			PrintSequencesClustal(rFile, alnLocal, uFlags);
			break;
		case FormatFasta:
		default:
			PrintSequencesFasta(rFile, alnLocal, uFlags);
	};
	if(g_csAfter != NULL)
		rFile << Subst(g_csAfter) << endl;
	return true;
}

// Fasta Output
void PrintSequencesFasta(ostream& os, const Tree::Alignment& aln, unsigned int)
{
	for(Tree::Alignment::const_iterator cit = aln.begin(); cit != aln.end(); ++cit)
	{
		// Print Label
		os << '>' << cit->first << endl;

		// Print Sequence in lines 60 charaters wide
		unsigned int u=0;
		while(u < cit->second.length())
		{
			os << cit->second[u];
			if((++u % 60) == 0)
				os << endl;
		}
		// Add newline if sequence did not fill up last row
		if(u%60)
			os << endl;
		os << endl;
	}
}

// Nexus Non-Interleaved Output
void PrintSequencesNexus(ostream &os, const Tree::Alignment& aln, unsigned int uFlags)
{
	// Write header for data block
	Tree::Alignment::const_iterator cit = aln.begin();
	os << "BEGIN DATA;" << endl;
	os << "\tDIMENSIONS NTAX=" << (unsigned int)aln.size() << " NCHAR=" << (unsigned int)cit->second.length() << ';' << endl;
	os << "\tFORMAT DATATYPE=";
	if(uFlags & FlagOutTranslate)
		os << "PROTEIN";
	else
		os << "DNA";
	os << " MISSING=? GAP=- MATCHCHAR=. EQUATE=\"+=- ==-\";" << endl;
	os << "\tMATRIX" << endl;
	
	// Write sequences in non-interleaved format
	for(; cit != aln.end(); ++cit)
		os << setw(15) << setiosflags(ios::left) << cit->first << ' ' << cit->second << endl;

	// Close data block
	os << ';' << endl << "END;" << endl << endl;
}

// Phylip Non-Interleaved Output
void PrintSequencesPhylip(ostream &os, const Tree::Alignment& aln, unsigned int)
{
	// Header
	Tree::Alignment::const_iterator cit = aln.begin();
	os << (unsigned int)aln.size() << ' ' << (unsigned int)cit->second.length() << endl;
	
	// Print non-interleaved sequences
	os << setfill(' ');
	for(; cit != aln.end(); ++cit)
		os << setw(10) << setiosflags(ios::left) << cit->first.substr(0,10) << cit->second << endl;
	os << endl;
}

// Clustal Output
void PrintSequencesClustal(ostream &os, const Tree::Alignment& aln, unsigned int)
{
	unsigned int uLen = (unsigned int)aln.begin()->second.length();
	unsigned int l;
	// Print interleaved sequences
	for(unsigned int u = 0; u < uLen; u+=l)
	{
		l = min(60u, uLen);
		// Print a row of each sequence
		for(Tree::Alignment::const_iterator cit = aln.begin(); cit != aln.end(); ++cit)
			os << setw(15) << setiosflags(ios::left) << cit->first << " " << cit->second.substr(u, l) << endl;
		os << endl << endl;
	}
}

// Replace extra gap characters with missing data
void FilterGapSingleChar(string& ss)
{
	int nGapState = 0;
	for(string::iterator it = ss.begin(); it!= ss.end(); ++it)
	{
		if(*it == '-' || *it == '=')
		{	
			if(nGapState == 1)
				*it = '?';
			else
				nGapState = 1; // Deletion
		}
		else if(*it == '+')
		{
			if(nGapState == 2)
				*it = '?';
			else
				nGapState = 2; // Insertion
		}
		else
			nGapState = 0;  // Bases
	}
}

// Replace special gap characters with '-'
void FilterGapNoPlus(string& ss)
{
	for(string::iterator it = ss.begin(); it!= ss.end(); ++it)
	{
		if(*it == '+' || *it == '=')
			*it = '-';
	}
}

// Convert uppercase sequences to lowercase sequences
void FilterLowerCase(string& ss)
{
	for(string::iterator it = ss.begin(); it!= ss.end(); ++it)
	{
		if(0x41 <= *it && *it <= 0x5A)
			*it |= 0x20;
	}
}

// DNA -> Protein
const char g_csAminoAcid[] = "KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG"; 

void FilterTranslate(string& ss)
{
	string::iterator it = ss.begin();
	string::const_iterator cit = ss.begin();
	while(*cit)
	{
		if(*cit < 'A')
		{
			*it++ = *cit;
			cit += 3;
		}
		else
		{
			// ACGT/acgt -> 0132
			// can be improved but requires endian-aware code

			unsigned int x = ((*cit++)&6u) << 3;
			x |= ((*cit++)&6u) << 1;
			x |= ((*cit++)&6u) >> 1;
			// Lookup letter from codon number
			*it++ = g_csAminoAcid[x];
		}
	}
	ss.resize(ss.size()/3);
}

// Filter the sequence based on the flags
void FilterSequence(const string& ssSrc, string& ssDest, unsigned int uFlags)
{
	// Copy Sequence
	ssDest = ssSrc;
	
	// Modify Sequence
	if(uFlags & FlagOutGapSingleChar)
		FilterGapSingleChar(ssDest);

	if(!(uFlags & FlagOutGapPlus))
		FilterGapNoPlus(ssDest);

	if(uFlags & FlagOutTranslate)
		FilterTranslate(ssDest);

	if(uFlags & FlagOutLowerCase)
		FilterLowerCase(ssDest);
}

// Filter the alignment based on the flags
void FilterAlignment(const Tree::Alignment& alnSrc, Tree::Alignment& alnDest, unsigned int uFlags)
{
	alnDest = alnSrc;
}
