// output.cc - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"

using namespace std;

void FilterSequence(const string& ssSrc, string& ssDest, unsigned long uFlags);
void PrintSequencesFasta(  ostream &os, const Tree::Alignment& aln, unsigned long uFlags);
void PrintSequencesNexus(  ostream &os, const Tree::Alignment& aln, unsigned long uFlags);
void PrintSequencesPhylip( ostream &os, const Tree::Alignment& aln, unsigned long uFlags);
void PrintSequencesClustal(ostream &os, const Tree::Alignment& aln, unsigned long uFlags);

FileFormat g_fileFormat = FASTA;
const char *g_csBlock = NULL;
int			g_nDataSet = 0;
int			g_nDataSetNum = 1;

// Reset Format
bool SetFormat(FileFormat fmt, int nNum, const char* csBlock)
{
	g_fileFormat = fmt;
	g_csBlock = csBlock;
	g_nDataSet = 0;
	g_nDataSetNum = nNum;
	return true;
}

// Initialize Output for formats that need it
void DawgIniOutput(ostream& os)
{
	switch(g_fileFormat)
	{
	case NEXUS:
		os << "#NEXUS" << endl << "[Created by DAWG Version " << VERSION << ']' << endl;
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

// Save alignment to output stream
bool SaveAlignment(ostream &rFile, const Tree::Alignment& aln, unsigned long uFlags)
{
	// Print DataSet number if multiple sequences will be returned
	if(g_nDataSetNum > 1)
		rFile << "[DataSet " << ++g_nDataSet << ']' << endl;

	// Filter sequences to a local alignment, applying format modifiers
	Tree::Alignment alnLocal;
	for(Tree::Alignment::const_iterator cit = aln.begin(); cit != aln.end(); ++cit)
		FilterSequence(cit->second, alnLocal[cit->first], uFlags);
	
	// Save Alignment in requested file format
	// Could be an array if flexibility is needed
	switch(g_fileFormat)
	{
		case NEXUS:
			PrintSequencesNexus(rFile, alnLocal, uFlags);
			break;
		case PHYLIP:
			PrintSequencesPhylip(rFile, alnLocal, uFlags);
			break;
		case CLUSTAL:
			PrintSequencesClustal(rFile, alnLocal, uFlags);
			break;
		case FASTA:
		default:
			PrintSequencesFasta(rFile, alnLocal, uFlags);
	};
	return true;
}

// Fasta Output
void PrintSequencesFasta(ostream& os, const Tree::Alignment& aln, unsigned long)
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
void PrintSequencesNexus(ostream &os, const Tree::Alignment& aln, unsigned long uFlags)
{
	// Write header for data block
	Tree::Alignment::const_iterator cit = aln.begin();
	os << "BEGIN DATA;" << endl;
	os << "\tDIMENSIONS NTAX=" << aln.size() << " NCHAR=" << cit->second.length() << ';' << endl;
	os << "\tFORMAT DATATYPE=";
	if(uFlags & FlagOutTranslate)
		os << "PROTEIN";
	else
		os << "DNA";
	os << " MISSING=? GAP=- MATCHCHAR=. EQUATE=\"+=-\";" << endl;
	os << "\tMATRIX" << endl;
	
	// Write sequences in non-interleaved format
	for(; cit != aln.end(); ++cit)
		os << setw(15) << setiosflags(ios::left) << cit->first << ' ' << cit->second << endl;

	// Close data block
	os << ';' << endl << "END;" << endl << endl;
	
	// Add Nexus Code Block if it exists
	if(g_csBlock)
		os << g_csBlock << endl << endl;
}

// Phylip Non-Interleaved Output
void PrintSequencesPhylip(ostream &os, const Tree::Alignment& aln, unsigned long)
{
	// Header
	Tree::Alignment::const_iterator cit = aln.begin();
	os << aln.size() << ' ' << cit->second.length() << endl;
	
	// Print non-interleaved sequences
	os << setfill(' ');
	for(; cit != aln.end(); ++cit)
		os << setw(10) << setiosflags(ios::left) << cit->first.substr(0,10) << cit->second << endl;
	os << endl;
}

// Clustal Output
void PrintSequencesClustal(ostream &os, const Tree::Alignment& aln, unsigned long)
{
	unsigned long uLen = aln.begin()->second.length();
	unsigned long l;
	// Print interleaved sequences
	for(unsigned long u = 0; u < uLen; u+=l)
	{
		l = min(60ul, uLen);
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
char g_csAminoAcid[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

void FilterTranslate(string& ss)
{
	// Setup a temporary string for the protein
	string ss2;
	ss2.reserve(ss.size()/3);
	
	// foreach codon
	for(string::iterator it = ss.begin(); it!= ss.end();)
	{
		// put a gap in the protein sequence if there is one in the nucleotide sequence
		if(*it == '-' || *it == '+' || *it == '=' || *it == '?')
		{
			ss2.push_back(*it);
			it += 3;
		}
		else
		{
			// determine the codon in [0,63]
			unsigned long uCodon = (Nucleotide::Letter2Number(*it++) << 4)
				| (Nucleotide::Letter2Number(*it++) << 2)
				| (Nucleotide::Letter2Number(*it++));
			// add the character that represents the codon to the temp string
			ss2.push_back( ((uCodon < 64) ? g_csAminoAcid[up] : '?'));
		}
	}
	// Set string to temp string
	ss = ss2;
}

// Filter the sequence based on the flags
void FilterSequence(const string& ssSrc, string& ssDest, unsigned long uFlags)
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
