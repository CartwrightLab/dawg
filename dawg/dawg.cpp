// dawg.cpp

#include "dawg.h"
#include "node.h"

#include <time.h>
#include <float.h>
#include <stdio.h>
#ifdef _WIN32
#	include <process.h>
#else
#	include <unistd.h>
#endif


#include <iostream>
#include <algorithm>
#include <iomanip>
#include <memory>

using namespace std;

bool Parse(const char* cs);
bool Execute();

// So we don't have to include "rand.h"
void mt_srand(unsigned long uKeys[], unsigned long uLen);
bool check_endian();

char csBuffer[32];

int main(int argc, char* argv[])
{
	if(!check_endian())
	{
		DawgError("Dawg compiled with the wrong endian settings.  See rand.h.");
		return 1;
	}
	bool bSerial = true; // Process files in serial fashion
	bool bUsage = (argc==1);
	bool bFullUsage = false;
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
					bFullUsage = true;
				case 'h':
				case 'H':
					bUsage = true;
					break;
				case 'u':
				case 'U':
					setvbuf(stdout, csBuffer, _IONBF, 32);
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
			cout << "Version " << DawgVersion() << endl << endl;
		}
		if(bUsage)
			cout << "Usage: dawg -[scvh?] file1 [file2 ...]" << endl;
		if(bFullUsage) { /*Todo*/ }

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
	int nTotalSeqLen = 0;
	double dGamma = 0.0, dIota = 0.0;
	double dNucFreq[4] = {0.25,0.25,0.25,0.25};
	double dRevParams[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	vector<double> vdParams;
	vector<double> vdRates;
	string ssSeq, ssModel = "JC", ssGapModel[2] = {"NB", "NB"};

	double dLambda[2] = {0.0, 0.0};
	vector<double> vdInsModel;
	vector<double> vdDelModel;

	vector<Tree*> vtTrees;
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
	
	if(!DawgVar::GetVector("Length", vSeqLen))
		vSeqLen.push_back(100);
    if(DawgVar::Get("Sequence", ssSeq))
	{
		if(vSeqLen.size() == 1)
			vSeqLen[0] = ssSeq.length();
		nTotalSeqLen = for_each(vSeqLen.begin(), vSeqLen.end(), SumValue<int>());
		if(nTotalSeqLen != (int)ssSeq.length())
			return DawgError("The sum of \"Length\" does not equal the length of \"Sequence\".");
	}
	if(vSeqLen.size() > 1 && vSeqLen.size() < vtTrees.size())
		return DawgError("When using variable section lengths, \"Length\" \
						 and \"Tree\" must have the same size.");
	DawgVar::GetVector("Rates", vdRates);  //bug (see below)
	if(vdRates.size() > 0 && vdRates.size() < (unsigned int)nTotalSeqLen)
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
	DawgVar* pVar = DawgVar::GetVar("GapParams");
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

	DawgVar::Get("File", ssFile);
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
	{
		mt_srand((unsigned long*)&vnSeed[0], vnSeed.size());
	}
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

	if(!Node::s_procSubst.Setup(dNucFreq, dRevParams))
		return DawgError("Invalid substitution model parameters.");
	
	if(!Nucleotide::Setup(dNucFreq, dGamma, dIota))
		return DawgError("Invalid G+I rates");
	
	//if(!Node::Scale(dScale))
	//	return DawgError("Invalid scaling parameter");
	Node::Scale(dScale);
	
	IndelModel::Params paramsDel, paramsIns;
	paramsIns.ssModel = ssGapModel[0];
	paramsIns.dLambda = dLambda[0];
	paramsIns.vdModel = vdInsModel;
	paramsDel.ssModel = ssGapModel[1];
	paramsDel.dLambda = dLambda[1];
	paramsDel.vdModel = vdDelModel;
	if(!Node::s_procIndel.Setup(paramsIns, paramsDel))
	
	
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
	
	SetFormat(fmt, nReps, ssBlock.empty() ? NULL : ssBlock.c_str(), bGapSingle);

	ofstream oFile;
	ostream* pOut;
	if(ssFile.empty())
		pOut = &cout;
	else if(DawgOpen(ssFile.c_str(), oFile))
		pOut = &oFile;
	else
		return DawgError("Unable to open \"%s\" for output.", ssFile.c_str());
	DawgIniOutput(*pOut);
	//Evolve
	while(nReps--)
	{
		vector<double>::iterator dit = vdRates.begin();
		for(unsigned int uTree = 0; uTree < vtTrees.size(); ++uTree)
		{
			Tree *pTree = vtTrees[uTree];
			Seq &rSeq = pTree->Sequence();
			int nSeqLen = (vSeqLen.size() == 1) ? vSeqLen[0]/vtTrees.size() : vSeqLen[uTree];

			rSeq.resize(nSeqLen);
			generate(rSeq.begin(), rSeq.end(), Nucleotide::Rand);
			if(!ssSeq.empty())
			{
				int nIndex = 0;
				for(Seq::iterator nit = rSeq.begin(); nit != rSeq.end(); ++nit)
					nit->m_nuc = CharToNuc(ssSeq[nIndex++]);
			}
			if(vdRates.size() > 0)
			{
				if(vdRates.size() < rSeq.size())
					return DawgError("Specified rate vector is smaller than sequence length.");
				for(Seq::iterator nit = rSeq.begin(); nit != rSeq.end(); ++nit)
					nit->m_dRate = *dit++;
			}
			pTree->ResetGaps();
			pTree->Evolve();
		}
		SaveSequences(*pOut, (const Node**)&(vtTrees[0]), vtTrees.size());
	}
	return true;
}
