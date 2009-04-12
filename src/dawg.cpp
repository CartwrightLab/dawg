/*  Dawg - DNA Assembly with Gaps - Simulating Sequence Evolution
    Copyright (c) 2004-2009  Reed A. Cartwright, PhD

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "dawg.h"
#include "tree.h"
#include "rand.h"
#include "var.h"

#include <getopt.h>

// Help Information
char g_csDawgTxt[] =
#include "dawgtxt.h"
;

// Notice
char g_csDawgNotice[] = 
#include "dawgnotice.h"
;

using namespace std;

bool Parse(const char* cs);
bool Execute();

bool g_bReportErrors = true;
bool g_bReportWarnings = true;

const char *g_csOutput = NULL;

int main(int argc, char* argv[])
{
	bool bSerial = true;

	// Check to see if Usage and Version need to be printed
	bool bUsage = (argc==1);
	bool bVersion = (argc==1);
	bool bOk = true;

	//Parse Aruments
	int ch;

	while((ch = getopt(argc, argv, "sScCvVhHuUbBqQeEwWo:O:")) != -1)
	{
		switch(ch)
		{
		//Serial Mode
		case 's':
		case 'S':
			bSerial = true;  
			break;
		//Combined Mode
		case 'c':
		case 'C':
			bSerial = false; 
			break;
		// Version Information
		case 'v':
		case 'V':
			bVersion = true;
			break;
		// Help Information
		case '?':
		case 'h':
		case 'H':
			bVersion = true;
			bUsage = true;
			break;
		// Unbuffered Output
		case 'u':
		case 'U':
			#ifdef SETVBUF_REVERSED
			setvbuf(stdout, _IONBF, NULL, 0);
			#else
			setvbuf(stdout, NULL, _IONBF, 0);
			#endif
			break;
		// Buffered Output
		case 'b':
		case 'B':
			#ifdef SETVBUF_REVERSED
			setvbuf(stdout, _IOFBF, NULL, BUFSIZ);
			#else
			setvbuf(stdout, NULL, _IOFBF, BUFSIZ);
			#endif
			break;
		// Disable Errors and Warnings
		case 'q':
		case 'Q':
			g_bReportErrors = false;
			g_bReportWarnings = false;
			break;
		// Enable Errors
		case 'e':
		case 'E':
			g_bReportErrors = true;
			break;
		// Enable Warnings
		case 'w':
		case 'W':
			g_bReportWarnings = true;
			break;
		// Output Override
		case 'o':
		case 'O':
			g_csOutput = optarg;
			break;
		// Error Reporting
		default:
			DawgError("Unreconized switch, \"%c\"", ch);
			bOk = false;
			bUsage = true;
			bVersion = true;
			break;
		}
	}
    argc -= optind;
    argv += optind;
	// Print Version and/or Usage Information and exit
	if(bVersion || bUsage)
	{
		if(bVersion)
		{
			cout << PACKAGE_STRING << endl 
				<< "DNA Assembly With Gaps" << endl
				<< "Copyright (C) 2004-2009 Reed A. Cartwright" << endl << endl
				<< g_csDawgNotice << endl
				<< "Send Bug Reports to " << PACKAGE_BUGREPORT << "." << endl << endl;
		}
		if(bUsage)
			cout << g_csDawgTxt << endl;

		return 0;
	}
	// Close if Error
	if(!bOk)
		return 1;
	if(bSerial)
	{
		// Process files in serial mode
		if(argc > 1)
			g_csOutput = NULL; // disable output override
		for(int i=0;i<argc;++i)
		{
			if(Parse(argv[i]))
			{
				if(!Execute())
					bOk = DawgError("Execution of \"%s\" failed.", argv[i]);
				DawgVar::ClearMap(); // Remove variables
			}
			else
				bOk = DawgError("Parsing of \"%s\" failed.", argv[i]);
		}
	}
	else
	{
		// Process files in combined mode
		for(int i=0;i<argc;++i)
		{
			if(!Parse(argv[i]))
				bOk = DawgError("Parsing of \"%s\" failed.", argv[i]);
		}
		if(!Execute())
			bOk = DawgError("Execution failed.");
	}
	return (int)!bOk;
}

// Fast and dirty random seed algorithm
// based off of NR's randqd1
inline unsigned int rand_seed()
{
	static unsigned int u = (unsigned int)(time(NULL)+3*getpid());
	return (u = u*1664525u + 1013904223u);
}

inline const char * SS2CS(const string& ss)
{
	return ss.empty() ? NULL : ss.c_str();
}

// Execute the Dawg process based on the current configuration
bool Execute()
{
	// Variables
	unsigned int uReps  =  1;
	vector<unsigned int> vuSeqLen;
	vector<string> vssSeqs;
	string::size_type uTotalSeqLen = 0, uTotalRateLen = 0;
	
	vector<double> vdGamma, vdIota, vdScale;

	double dNucFreq[4] = {0.25,0.25,0.25,0.25};
	double dRevParams[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	vector<double> vdParams;
	vector< vector<double> > vvdRates;
	string ssModel = "JC", ssGapModel[2] = {"US", "US"};

	double dLambda[2] = {0.0, 0.0};
	vector<double> vdGapModel[2];
	vdGapModel[0].push_back(1.0);
	vdGapModel[1].push_back(1.0);
	unsigned int uKeepFlank = 0;

	vector<NewickNode*> vtTrees;
	double	dTreeScale = 1.0;
	vector<unsigned int>   vuSeed;

	string ssFile = "-";
	unsigned int uFmt = FormatClustal;
	string ssFormat = "";
	string ssOutBlockHead = "";
	string ssOutBlockBefore = "";
	string ssOutBlockAfter = "";
	string ssOutBlockTail = "";
	bool bOutSubst = true;

	bool bGapSingle = false, bGapPlus = false, bLowerCase = false, bTranslate = false;
	bool bKeepEmpty = false;
	unsigned int uWidth = 1;
	DawgVar::Vec::size_type nRes;

	// Read variables from configuration

	if(!DawgVar::GetVector("Tree", vtTrees))
		return DawgError("No trees specified.");

	if(DawgVar::GetVector("Sequence", vssSeqs))
	{
		if(vssSeqs.size() < vtTrees.size())
			return DawgError("\"Sequence\" and \"Tree\" must have the same size.");
		for(vector<string>::const_iterator cit = vssSeqs.begin(); cit != vssSeqs.end(); ++cit)
			uTotalSeqLen += cit->length();
	}
	else if(DawgVar::GetVector("Length", vuSeqLen))
	{
		if(vuSeqLen.size() < vtTrees.size())
			return DawgError("\"Length\" and \"Tree\" must have the same size.");
		for(vector<unsigned int>::const_iterator cit = vuSeqLen.begin(); cit != vuSeqLen.end(); ++cit)
			uTotalSeqLen += *cit;
	}
	else
		vuSeqLen.resize(vtTrees.size(), 100);
	
	DawgVar::Get("Width", uWidth);

	vdGamma.resize(uWidth, 0.0);
	vdIota.resize(uWidth, 0.0);
	vdScale.resize(uWidth, 1.0);

	vvdRates.resize(vtTrees.size());
	if(DawgVar::GetMatrix("Rates", &vvdRates[0], vvdRates.size()))
	{
		for(vector< vector<double> >::const_iterator cit = vvdRates.begin();
			cit != vvdRates.end(); ++cit)
			uTotalRateLen += cit->size();
		if(uTotalRateLen > 0 && uTotalRateLen < uTotalSeqLen)
			return DawgError("\"Rates\" vector is too small");
	}
	else
		vvdRates.clear();

    DawgVar::Get("Reps", uReps);
	DawgVar::GetVector("Seed", vuSeed);
	DawgVar::Get("Model", ssModel);
	DawgVar::GetVector("Params", vdParams);
	if(!DawgVar::GetArray("Gamma", &vdGamma[0], uWidth)) // Coef of Variation
	{
		if(DawgVar::GetArray("Alpha", &vdGamma[0], uWidth))  // Shape parameter
		{
			for(vector<double>::iterator it = vdGamma.begin();
				it != vdGamma.end(); ++it)
			{
				*it = 1.0 / *it;
			}
		}
	}
	DawgVar::GetArray("Iota", &vdIota[0], uWidth);
	DawgVar::GetArray("Scale", &vdScale[0], uWidth);
	DawgVar::Get("TreeScale", dTreeScale);
	DawgVar::GetArray("GapModel", ssGapModel, 2);
	nRes = DawgVar::GetArray("Freqs", dNucFreq, 4, false);
	if(nRes > 0 && nRes < 4)
		return DawgError("\"Freqs\" specified incorrectly.");
	
	DawgVar::Get("GapSingleChar", bGapSingle);
	DawgVar::Get("GapPlus", bGapPlus);
	DawgVar::Get("LowerCase", bLowerCase);
	DawgVar::Get("Translate", bTranslate);
	DawgVar::Get("KeepEmpty", bKeepEmpty);

	nRes = DawgVar::GetArray("Lambda", dLambda, 2, false);
	if(nRes)
	{
		if(nRes == 1)
			dLambda[1] = dLambda[0] *= 0.5;
		nRes = DawgVar::GetMatrix("GapParams", vdGapModel, 2);
		if(nRes == 0)
			return DawgError("\"GapParams\" must be specified if \"Lambda\" is.");
	}
	DawgVar::Get("KeepFlank", uKeepFlank);

	DawgVar::Get("File", ssFile);
	DawgVar::Get("Format", ssFormat);
    if(DawgVar::Get("NexusCode", ssOutBlockAfter))
		DawgWarn("NexusCode is depreciated.  Use Out.Block.* instead.");
	
	DawgVar::Get("Out.Block.Head", ssOutBlockHead);
	DawgVar::Get("Out.Block.Before", ssOutBlockBefore);
	DawgVar::Get("Out.Block.After", ssOutBlockAfter);
	DawgVar::Get("Out.Block.Tail", ssOutBlockTail);
	DawgVar::Get("Out.Subst", bOutSubst);

	// Setup Model Parameters from Variables

	// Load random number generator
	if(vuSeed.empty())
	{
		for(unsigned int u=0;u<4;++u)
			vuSeed.push_back(rand_seed());
	}
	mt_srand(&vuSeed[0], vuSeed.size());
	
	// Setup recombinant tree
	Tree myTree;
	myTree.ProcessTree(vtTrees.begin(), vtTrees.end());

	// Construct substitution matrices
	if(ssModel == "GTR")
	{
		if(vdParams.size() < 6)
			return DawgError("GTR model requires at least 6 numerical parameters.");
		for(int i=0;i<6;i++)
			dRevParams[i] = vdParams[i];
	}
	else if(ssModel == "JC")
	{
		dNucFreq[3] = dNucFreq[2] = dNucFreq[1] = dNucFreq[0] = 0.25;
		dRevParams[5] = dRevParams[4] = dRevParams[3] = 1.0;
		dRevParams[2] = dRevParams[1] = dRevParams[0] = 1.0;
	}
	else if(ssModel == "K2P")
	{
		if(vdParams.size() < 2)
			return DawgError("K2P model requires at least 2 numerical parameters.");
		dNucFreq[3] = dNucFreq[2] = dNucFreq[1] = dNucFreq[0] = 0.25;
		dRevParams[4] = dRevParams[1] = vdParams[0];
		dRevParams[0] = dRevParams[2] = dRevParams[3] = dRevParams[5] = vdParams[1];
	}
	else if(ssModel == "K3P")
	{
		if(vdParams.size() < 3)
			return DawgError("K3P model requires at least 3 numerical parameters.");
		dNucFreq[3] = dNucFreq[2] = dNucFreq[1] = dNucFreq[0] = 0.25;
		dRevParams[4] = dRevParams[1] = vdParams[0];
		dRevParams[2] = dRevParams[3] = vdParams[1];
		dRevParams[0] = dRevParams[5] = vdParams[2];

	}
	else if(ssModel == "HKY")
	{
		if(vdParams.size() < 1)
			return DawgError("HKY model requires at least 1 numerical parameter.");
		dRevParams[4] = dRevParams[1] = vdParams[0];
		dRevParams[0] = dRevParams[2] = dRevParams[3] = dRevParams[5] = 1.0;
	}
	else if(ssModel == "F81")
	{
		dRevParams[5] = dRevParams[4] = dRevParams[3] = 1.0;
		dRevParams[2] = dRevParams[1] = dRevParams[0] = 1.0;
	}
	else if(ssModel == "F84")
	{
		if(vdParams.size() < 1)
			return DawgError("F84 model requires at least 1 numerical parameter.");
		dRevParams[1] = 1.0 + vdParams[0]/(dNucFreq[0]+dNucFreq[2]);
		dRevParams[4] = 1.0 + vdParams[0]/(dNucFreq[1]+dNucFreq[3]);
		dRevParams[0] = dRevParams[2] = dRevParams[3] = dRevParams[5] = 1.0;
	}
	else if(ssModel == "TN")
	{
		if(vdParams.size() < 3)
			return DawgError("TN model requires at least 3 numerical parameters.");
		dRevParams[1] = vdParams[0];
		dRevParams[4] = vdParams[1];
		dRevParams[0] = dRevParams[2] = dRevParams[3] = dRevParams[5] = vdParams[2];
	}
	else
		return DawgError("Unknown substitution model, \"%s\".", ssModel.c_str());
	
	// Construct Indel parameters
	IndelModel::Params paramsDel, paramsIns;
	paramsIns.ssModel = ssGapModel[0];
	paramsIns.dLambda = dLambda[0];
	paramsIns.vdModel = vdGapModel[0];

	paramsDel.ssModel = ssGapModel[1];
	paramsDel.dLambda = dLambda[1];
	paramsDel.vdModel = vdGapModel[1];
	
	// Initialize Evolution
	if(!myTree.SetupEvolution(dNucFreq, dRevParams, paramsIns, paramsDel,
		uWidth, vdGamma, vdIota, vdScale, dTreeScale, uKeepFlank ))
		return DawgError("Bad evolution parameters");

	// Initialize Root
	if(!myTree.SetupRoot(vssSeqs, vuSeqLen, vvdRates))
		return DawgError("Bad root parameters");
	
	// Read Out.Block.* from file if neccessary
   	if(!ssOutBlockHead.empty())
	{
		ifstream iFile(ssOutBlockHead.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockHead, '\0');
	}
   	if(!ssOutBlockBefore.empty())
	{
		ifstream iFile(ssOutBlockBefore.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockBefore, '\0');
	}
   	if(!ssOutBlockAfter.empty())
	{
		ifstream iFile(ssOutBlockAfter.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockAfter, '\0');
	}
   	if(!ssOutBlockTail.empty())
	{
		ifstream iFile(ssOutBlockTail.c_str());
		if(iFile.is_open())
			getline(iFile, ssOutBlockTail, '\0');
	}

	// Output Format
	if(ssFormat == "")
	{
		string::size_type pos = ssFile.find_first_of(':');
		if( pos != string::npos)
			ssFormat = ssFile.substr(0, pos);
		else
		{
			pos = ssFile.find_last_of('.');
			if(pos != string::npos)
				ssFormat = ssFile.substr(pos+1);
			else
				ssFormat = "Clustal";
		}
	}
	// Override output
	if(g_csOutput != NULL)
	{
		ssFile = g_csOutput;
		string::size_type pos = ssFile.find_first_of(':');
		if( pos != string::npos)
			ssFormat = ssFile.substr(0, pos);
		else
		{
			pos = ssFile.find_last_of('.');
			if(pos != string::npos)
				ssFormat = ssFile.substr(pos+1);
		}
	}
	if(ssFormat == "Fasta" || ssFormat == "fas")
		uFmt = FormatFasta;
	else if(ssFormat == "Nexus" || ssFormat == "nex")
		uFmt = FormatNexus;
	else if(ssFormat == "Phylip" || ssFormat == "phy")
		uFmt = FormatPhylip;
	else if(ssFormat == "Clustal" || ssFormat == "aln"
		|| ssFormat == "poo" || ssFormat == "txt" || ssFormat == "out")
		uFmt = FormatClustal;
	else
		return DawgError("Unknown file format, \"%s\".", ssFormat.c_str());

	SetFormat(uFmt, uReps, SS2CS(ssOutBlockHead), SS2CS(ssOutBlockBefore),
		SS2CS(ssOutBlockAfter), SS2CS(ssOutBlockTail), bOutSubst);
	
	// Check translate parameter
	if(bTranslate && uWidth != 3)
		return DawgError("Translate requires a Width of 3.");
	
	// setup output flags
	unsigned int uOutFlags = 0u;
	if(bGapSingle)
		uOutFlags |= FlagOutGapSingleChar;
	if(bGapPlus)
		uOutFlags |= FlagOutGapPlus|FlagOutKeepEmpty;
	if(bLowerCase)
		uOutFlags |= FlagOutLowerCase;
	if(bTranslate)
		uOutFlags |= FlagOutTranslate;
	if(bKeepEmpty)
		uOutFlags |= FlagOutKeepEmpty;

	// setup output location
	ostream* pOut;
	ofstream ofOut;
	string::size_type pos = ssFile.find_first_of(':');
	if( pos != string::npos)
		ssFile.erase(0, pos+1);

	if(ssFile == "-" || ssFile.empty())
		pOut = &cout;
	else
	{
		ofOut.open(ssFile.c_str());
		if(!ofOut.is_open())
			return DawgError("Unable to open \"%s\" for output.", ssFile.c_str());
		pOut = &ofOut;
	}
	DawgIniOutput(*pOut);

	// Evolve Many Sequences
	while(uReps--)
	{
		//Evolve
		myTree.Evolve();

		//SaveOutput
		Tree::Alignment aln;
		myTree.Align(aln, uOutFlags);
		if(!SaveAlignment(*pOut, aln, uOutFlags))
			return DawgError("Error saving alignment.");
	}

	DawgFinOutput(*pOut);

	return true;
}

// Error Reporting Routine
bool DawgError(const char* csErr, ...)
{
	if(!g_bReportErrors)
		return false;
	fprintf(stderr, "Error: ");
	va_list args;
	va_start(args, csErr);
#if defined(HAVE_VPRINTF)
	vfprintf(stderr, csErr, args);
#elif defined(HAVE_DOPRNT)
	_doprnt(csErr, args, stderr);
#else
	fprintf(stderr, "%s", csErr);
#endif
	va_end(args);
	fprintf(stderr, "\n");
	return false;
}

bool DawgWarn(const char* csErr, ...)
{
	if(!g_bReportWarnings)
		return false;
	fprintf(stderr, "Warning: ");
	va_list args;
	va_start(args, csErr);
#if defined(HAVE_VPRINTF)
	vfprintf(stderr, csErr, args);
#elif defined(HAVE_DOPRNT)
	_doprnt(csErr, args, stderr);
#else
	fprintf(stderr, "%s", csErr);
#endif
	va_end(args);
	fprintf(stderr, "\n");
	return false;
}
