// dawg.cc - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

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
			cout << "Copyright (C) 2004 Reed A. Cartwright (all rights reserved)" << endl;
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
	unsigned long uReps  =  1;
	vector<unsigned long> vuSeqLen;
	vector<string> vssSeqs;
	unsigned long uTotalSeqLen = 0, uTotalRateLen = 0;
	
	vector<double> vdGamma, vdIota, vdScale;

	double dNucFreq[4] = {0.25,0.25,0.25,0.25};
	double dRevParams[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	vector<double> vdParams;
	vector< vector<double> > vvdRates;
	string ssModel = "JC", ssGapModel[2] = {"NB", "NB"};

	double dLambda[2] = {0.0, 0.0};
	vector<double> vdGapModel[2];

	vector<NewickNode*> vtTrees;
	double	dTreeScale = 1.0;
	vector<unsigned long>   vuSeed;



	string ssFile = "-";
	FileFormat fmt = FASTA;
	string ssFormat;
	string ssNexusCode;

	bool bGapSingle = false, bGapPlus = false;
	unsigned long uWidth = 1;
	int nRes;

	// Ready Variables

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
		for(vector<unsigned long>::const_iterator cit = vuSeqLen.begin(); cit != vuSeqLen.end(); ++cit)
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

	nRes = DawgVar::GetArray("Lambda", dLambda, 2);
	if(nRes)
	{
		nRes = DawgVar::GetMatrix("GapParams", vdGapModel, 2);
		if(nRes == 0)
			return DawgError("\"GapParams\" must be specified if \"Lambda\" is.");
	}

	DawgVar::Get("File", ssFile);
	DawgVar::Get("Format", ssFormat);
	DawgVar::Get("NexusCode", ssNexusCode);

	// Load Variables
	if(vuSeed.empty())
	{
		for(unsigned int u=0;u<4;++u)
			vuSeed.push_back(rand_seed());
	}
	mt_srand(&vuSeed[0], vuSeed.size());
	
	Tree myTree;
	for(vector<NewickNode*>::const_iterator treeit = vtTrees.begin(); treeit != vtTrees.end(); ++treeit)
		myTree.ProcessTree(*treeit);

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
	paramsIns.vdModel = vdGapModel[0];

	paramsDel.ssModel = ssGapModel[1];
	paramsDel.dLambda = dLambda[1];
	paramsDel.vdModel = vdGapModel[1];
	
	if(!myTree.SetupEvolution(dNucFreq, dRevParams, paramsIns, paramsDel,
		uWidth, vdGamma, vdIota, vdScale, dTreeScale ))
		return DawgError("Bad evolution parameters");
	if(!myTree.SetupRoot(vssSeqs, vuSeqLen, vvdRates))
		return DawgError("Bad root parameters");
	
   	if(!ssNexusCode.empty())
	{
		ifstream iFile(ssNexusCode.c_str());
		if(iFile.is_open())
			getline(iFile, ssNexusCode, '\0');
	}
	if(ssFormat.empty() || ssFormat == "Fasta")
		fmt = FASTA;
	else if(ssFormat == "Nexus")
		fmt = NEXUS;
	else if(ssFormat == "Phylip")
		fmt = PHYLIP;
	else
		return DawgError("Unknown file format, \"%s\".");
	
	SetFormat(fmt, uReps, ssNexusCode.empty() ? NULL : ssNexusCode.c_str());
	
	ostream* pOut;
	ofstream ofOut;
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

	while(uReps--)
	{
		//Evolve
		myTree.Evolve();

		//SaveOutput
		Tree::Alignment aln;
		myTree.Align(aln, bGapPlus, bGapSingle);
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
