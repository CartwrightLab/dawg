// dawg.cpp

#include "dawg.h"
#include "tree.h"
#include "rand.h"
#include "var.h"

char g_csDawgTxt[] =
#include "dawgtxt.h"
;

using namespace std;

bool Parse(const char* cs);
bool Execute();

char csBuffer[32];

int main(int argc, char* argv[])
{
	bool bSerial = true; // Process files in serial fashion
	bool bUsage = (argc==1);
	bool bVersion = (argc==1);
	bool bOk = true;
	//Parse Cmds
	int i=1;
	for(;i<argc;++i)
	{
		char* pch = argv[i];
		if( *pch != '-' || *(pch+1) == '\0')
			break;
		while(*++pch)
		{
			switch(*pch)
			{
				case 's':
				case 'S':
					bSerial = true;  //Serial Mode
					break;
				case 'c':
				case 'C':
					bSerial = false; //Combine Mode
					break;
				case 'v':
				case 'V':
					bVersion = true;
					break;
				case '?':
				case 'h':
				case 'H':
					bVersion = true;
					bUsage = true;
					break;
				case 'u':
				case 'U':
					#ifdef SETVBUF_REVERSED
					setvbuf(stdout, _IONBF, csBuffer, 32);
					#else
					setvbuf(stdout, csBuffer, _IONBF, 32);
					#endif
					break;
				default:
					DawgError("Unreconized switch, \"%c\"", *pch);
					bOk = false;
					bUsage = true;
					bVersion = true;
					break;
			};
		}
	}
	if(bVersion || bUsage)
	{
		if(bVersion)
		{
			cout << "DAWG - DNA Assembly With Gaps" << endl;
			cout << "Copyright (c) 2003-2004 Reed A. Cartwright (all rights reserved)" << endl;
			cout << "Version " << VERSION << endl << endl;
		}
		if(bUsage)
			cout << g_csDawgTxt << endl;

		return 0;
	}
	if(!bOk)
		return 1;
	if(bSerial)
	{
		for(;i<argc;++i)
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
		for(;i<argc;++i)
		{
			if(!Parse(argv[i]))
				bOk = DawgError("Parsing of \"%s\" failed.", argv[i]);
		}
		if(!Execute())
			bOk = DawgError("Execution failed.");
	}
	return (int)!bOk;
}

inline unsigned int rand_seed()
{
	// based off of NR's randqd1	
	static unsigned long u = time(NULL)+3*getpid();
	return (u = u*1664525u + 1013904223u);
}

bool Execute()
{
	// Variables
	int nReps  =  1;
	vector<int> vSeqLen;
	vector<string> vSeqs;
	int nTotalSeqLen = 0, nTotalRateLen = 0;
	double dGamma = 0.0, dIota = 0.0;
	double dNucFreq[4] = {0.25,0.25,0.25,0.25};
	double dRevParams[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	vector<double> vdParams;
	vector< vector<double> > vvdRates;
	string ssModel = "JC", ssGapModel[2] = {"NB", "NB"};

	double dLambda[2] = {0.0, 0.0};
	vector<double> vdInsModel;
	vector<double> vdDelModel;

	vector<NewickNode*> vtTrees;
	double		  dScale = 1.0;
	vector<int>   vnSeed;

	string ssFile;
	FileFormat fmt = FASTA;
	string ssFormat;
	string ssBlockFile;
	string ssBlock;

	bool bGapSingle = false;

	// Ready Variables

	if(!DawgVar::GetVector("Tree", vtTrees))
		return DawgError("No trees specified.");

	if(DawgVar::GetVector("Sequence", vSeqs))
	{
		if(vSeqs.size() < vtTrees.size())
			return DawgError("\"Sequence\" and \"Tree\" must have the same size.");
		for(vector<string>::const_iterator cit = vSeqs.begin(); cit != vSeqs.end(); ++cit)
			nTotalSeqLen += cit->length();
	}
	else if(DawgVar::GetVector("Length", vSeqLen))
	{
		if(vSeqLen.size() < vtTrees.size())
			return DawgError("\"Length\" and \"Tree\" must have the same size.");
		for(vector<int>::const_iterator cit = vSeqLen.begin(); cit != vSeqLen.end(); ++cit)
			nTotalSeqLen += *cit;
	}
	else
		vSeqLen.resize(vtTrees.size(), 100);

	DawgVar *pVar = DawgVar::GetVar("Rates");
	if(pVar != NULL && pVar->Size())
	{
		double dTemp;
		vvdRates.resize(pVar->Size());
		if(pVar->IsType(DawgVar::tyVector))
		{
			if(pVar->GetAt(0).IsType(DawgVar::tyVector))
			{
				// Assume Matrix
				for(unsigned int u=0; u<pVar->Size(); ++u)
					pVar->GetAt(u).GetVector(vvdRates[u]);
			}
			else
			{
				// Assume Vector
				if(!pVar->GetVector(vvdRates[0]))
					return DawgError("Error reading \"Rates\" vector.");
			}
		}
		else if(pVar->Get(dTemp))
		{
			vvdRates[0].push_back(dTemp);
		}
		else
			return DawgError("Error reading \"Rates\" vector.");
		for(vector< vector<double> >::const_iterator cit = vvdRates.begin();
			cit != vvdRates.end(); ++cit)
			nTotalRateLen += cit->size();
	}
	if(vvdRates.size() > 0 && nTotalRateLen < nTotalSeqLen)
		return DawgError("\"Rates\" vector is too small");

    DawgVar::Get("Reps", nReps);
	DawgVar::GetVector("Seed", vnSeed);
	DawgVar::Get("Model", ssModel);
	DawgVar::GetVector("Params", vdParams);
	if(!DawgVar::Get("Gamma", dGamma)) // Coef of Variation
	{
		if(DawgVar::Get("Alpha", dGamma))  // Shape parameter
			dGamma = 1.0/dGamma;
	}
	DawgVar::Get("Iota", dIota);
	DawgVar::Get("Scale", dScale); // make a vector (?)
	int nRes = DawgVar::GetArray("GapModel", ssGapModel, 2);
	if(nRes >0 && nRes < 2)
		ssGapModel[2] = ssGapModel[1];

	nRes = DawgVar::GetArray("Freqs", dNucFreq, 4);
	if(nRes > 0 && nRes < 4)
		return DawgError("\"Freqs\" specified incorrectly.");
	
	DawgVar::Get("GapSingleChar", bGapSingle);

	nRes = DawgVar::GetArray("Lambda", dLambda, 2);
	if(nRes == 1)
		dLambda[1] = dLambda[0];
	pVar = DawgVar::GetVar("GapParams");
	if(pVar == NULL || pVar->Size() == 0)
	{
		if(nRes)
			return DawgError("\"GapParams\" must be specified if \"Lambda\" is.");
	}
	else if((*pVar)[0].IsType(DawgVar::tyVector))
	{
		(*pVar)[0].GetVector(vdInsModel);
		if(pVar->Size() == 1)
			vdDelModel = vdInsModel;
		else
			(*pVar)[1].GetVector(vdDelModel);
	}
	else
	{
		pVar->GetVector(vdInsModel);
		vdDelModel = vdInsModel;
	}

	if(!DawgVar::Get("File", ssFile))
		ssFile = "-";
	DawgVar::Get("Format", ssFormat);
	if(!DawgVar::Get("NexusBlock", ssBlock))
		DawgVar::Get("NexusBlockFile", ssBlockFile);

	// Load Variables
	if(vnSeed.empty())
	{
		unsigned long uSeed[4];
		uSeed[0] = rand_seed();
		uSeed[1] = rand_seed();
		uSeed[2] = rand_seed();
		uSeed[3] = rand_seed();
		mt_srand(uSeed, 4);
	}
	else
		mt_srand((unsigned long*)&vnSeed[0], vnSeed.size());
	
	Tree myTree;
	for(vector<NewickNode*>::const_iterator treeit = vtTrees.begin(); treeit != vtTrees.end(); ++treeit)
		myTree.ProcessTree(*treeit);

	if(ssModel == "GTR")
	{
		if(vdParams.size() < 6)
			return DawgError("REV model requires at least 6 numerical parameters.");
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
		return DawgError("Unknown Model, \"%s\"", ssModel.c_str());
	
	IndelModel::Params paramsDel, paramsIns;
	paramsIns.ssModel = ssGapModel[0];
	paramsIns.dLambda = dLambda[0];
	paramsIns.vdModel = vdInsModel;
	paramsDel.ssModel = ssGapModel[1];
	paramsDel.dLambda = dLambda[1];
	paramsDel.vdModel = vdDelModel;
	
	if(!myTree.SetupEvolution(dNucFreq, dRevParams, paramsIns, paramsDel,
		dGamma, dIota, dScale))
		return DawgError("Bad evolution parameters");
	if(!myTree.SetupRoot(vSeqs, vSeqLen, vvdRates))
		return DawgError("Bad root parameters");
	
   	if(ssBlock.empty() && !ssBlockFile.empty())
	{
		ifstream iFile(ssBlockFile.c_str());
		if(!iFile.is_open())return DawgError("Unable to open \"%s\" for block file.", ssBlockFile.c_str());
		getline(iFile, ssBlock, '\0');
	}
	if(ssFormat.empty() || ssFormat == "Fasta")
		fmt = FASTA;
	else if(ssFormat == "Nexus")
		fmt = NEXUS;
	else if(ssFormat == "Phylip")
		fmt = PHYLIP;
	else
		return DawgError("Unknown file format, \"%s\".");
	
	SetFormat(fmt, nReps, ssBlock.empty() ? NULL : ssBlock.c_str());
	
	ostream* pOut;
	ofstream ofOut;
	if(ssFile.empty() || ssFile == "-")
		pOut = &cout;
	else
	{
		ofOut.open(ssFile.c_str());
		if(!ofOut.is_open())
			return DawgError("Unable to open \"%s\" for output.", ssFile.c_str());
		pOut = &ofOut;
	}
	DawgIniOutput(*pOut);

	while(nReps--)
	{
		//Evolve
		myTree.Evolve();

		//SaveOutput
		Tree::Alignment aln;
		myTree.Align(aln, bGapSingle);
		if(!SaveAlignment(*pOut, aln))
			return DawgError("Error saving alignment.");
	}

	return true;
}

bool DawgError(const char* csErr, ...)
{
	fprintf(stderr, "Error: ");
	va_list args;
	va_start(args, csErr);
	vfprintf(stderr, csErr, args);
	va_end(args);
	fprintf(stderr, "\n");
	return false;
}
