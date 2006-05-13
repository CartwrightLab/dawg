#ifndef DAWG_DAWGVAR_H
#define DAWG_DAWGVAR_H

#include <string>
#include <vector>
#include <fstream>
#include "var.h"

namespace Dawg {

class Variables
{
public:
	std::vector< std::vector< double > > dIndParams;
	std::vector< unsigned int > uSeqSeed;
	std::vector< double > dSeqTreeScale;
	bool bOutSeqGapPlus;
	std::string ssOutBlockHead;
	std::string ssOutFormat;
	std::string ssIndDelModel;
	bool bOutSeqTranslate;
	bool bOutSeqLowerCase;
	bool bOutReplace;
	std::vector< double > dIndDelParams;
	double dSubFreqs[4];
	std::vector< double > dIndInsParams;
	double dIndDelLambda;
	std::vector< int > nSeqLength;
	std::string ssIndInsModel;
	double dIndLambda[2];
	double dIndInsLambda;
	bool bOutSeqGapSingleChar;
	std::vector< double > dSubIota;
	std::vector< double > dSubAlpha;
	double dSeqReps;
	std::string ssOutBlockTail;
	std::string ssOutFile;
	std::vector< double > dSubScale;
	std::vector< NewickNode * > pSeqTree;
	int nSubCodonWidth;
	std::vector< double > dSubGamma;
	std::string ssIndModel[2];
	std::vector< std::vector< double > > dSeqRates;
	std::string ssOutBlockAfter;
	double dSubParams[6];
	std::string ssSubModel;
	std::vector< std::string > ssSeqSequence;
	std::string ssOutBlockBefore;

	bool ParseFile(const char* cs)
	{
		ToDefault();
		VarDB db;
		if(!db.Parse(cs))
			return false;
		(db.Get("Ind.Params", dIndParams)
			|| db.Get("GapParams", dIndParams));
		(db.Get("Seq.Seed", uSeqSeed)
			|| db.Get("Seed", uSeqSeed));
		(db.Get("Seq.TreeScale", dSeqTreeScale)
			|| db.Get("TreeScale", dSeqTreeScale));
		(db.Get("Out.Seq.GapPlus", bOutSeqGapPlus)
			|| db.Get("GapPlus", bOutSeqGapPlus));
		(db.Get("Out.Block.Head", ssOutBlockHead));
		LoadStringOrFile(ssOutBlockHead);
		(db.Get("Out.Format", ssOutFormat));
		(db.Get("Ind.Del.Model", ssIndDelModel));
		(db.Get("Out.Seq.Translate", bOutSeqTranslate)
			|| db.Get("Translate", bOutSeqTranslate));
		(db.Get("Out.Seq.LowerCase", bOutSeqLowerCase)
			|| db.Get("LowerCase", bOutSeqLowerCase));
		(db.Get("Out.Replace", bOutReplace)
			|| db.Get("Out.Subst", bOutReplace));
		(db.Get("Ind.Del.Params", dIndDelParams));
		(db.Get("Sub.Freqs", dSubFreqs)
			|| db.Get("Freqs", dSubFreqs));
		(db.Get("Ind.Ins.Params", dIndInsParams));
		(db.Get("Ind.Del.Lambda", dIndDelLambda));
		(db.Get("Seq.Length", nSeqLength)
			|| db.Get("Length", nSeqLength));
		(db.Get("Ind.Ins.Model", ssIndInsModel));
		(db.Get("Ind.Lambda", dIndLambda)
			|| db.Get("Lambda", dIndLambda));
		(db.Get("Ind.Ins.Lambda", dIndInsLambda));
		(db.Get("Out.Seq.GapSingleChar", bOutSeqGapSingleChar)
			|| db.Get("GapSingleChar", bOutSeqGapSingleChar));
		(db.Get("Sub.Iota", dSubIota)
			|| db.Get("Iota", dSubIota));
		(db.Get("Sub.Alpha", dSubAlpha)
			|| db.Get("Alpha", dSubAlpha));
		(db.Get("Seq.Reps", dSeqReps)
			|| db.Get("Reps", dSeqReps));
		(db.Get("Out.Block.Tail", ssOutBlockTail));
		LoadStringOrFile(ssOutBlockTail);
		(db.Get("Out.File", ssOutFile));
		(db.Get("Sub.Scale", dSubScale)
			|| db.Get("Scale", dSubScale));
		if(!(db.Get("Seq.Tree", pSeqTree)
			|| db.Get("Tree", pSeqTree)))
			return DawgError("Seq.Tree is not set properly.  It needs to be a 'Vector of Trees'.");
		(db.Get("Sub.CodonWidth", nSubCodonWidth)
			|| db.Get("Nuc.Width", nSubCodonWidth)
			|| db.Get("Width", nSubCodonWidth));
		(db.Get("Sub.Gamma", dSubGamma)
			|| db.Get("Gamma", dSubGamma));
		(db.Get("Ind.Model", ssIndModel)
			|| db.Get("GapModel", ssIndModel));
		(db.Get("Seq.Rates", dSeqRates)
			|| db.Get("Rates", dSeqRates));
		(db.Get("Out.Block.After", ssOutBlockAfter)
			|| db.Get("NexusCode", ssOutBlockAfter));
		LoadStringOrFile(ssOutBlockAfter);
		(db.Get("Sub.Params", dSubParams)
			|| db.Get("Params", dSubParams));
		(db.Get("Sub.Model", ssSubModel)
			|| db.Get("Model", ssSubModel));
		(db.Get("Seq.Sequence", ssSeqSequence)
			|| db.Get("Sequence", ssSeqSequence));
		(db.Get("Out.Block.Before", ssOutBlockBefore));
		LoadStringOrFile(ssOutBlockBefore);
		
		return true;
	}
	void ToDefault()
	{
		static double dIndParamsDefault[2] = {1.5, 1000};
		dIndParams = std::vector< std::vector< double > >(2, std::vector<double>(dIndParamsDefault, dIndParamsDefault+2));
		uSeqSeed = std::vector< unsigned int >(1,0);
		dSeqTreeScale = std::vector< double >(1,1.0);
		bOutSeqGapPlus = false;
		ssOutBlockHead = "";
		ssOutFormat = "Clustal";
		ssIndDelModel = "PL";
		bOutSeqTranslate = false;
		bOutSeqLowerCase = false;
		bOutReplace = false;
		dIndDelParams = std::vector< double >(dIndParamsDefault, dIndParamsDefault+2);
		static double dSubFreqsDefault[4] = {0.25, 0.25, 0.25, 0.25};
		memcpy(dSubFreqs, dSubFreqsDefault, sizeof(dSubFreqs));
		dIndInsParams = std::vector< double >(dIndParamsDefault, dIndParamsDefault+2);
		dIndDelLambda = 0.0;
		nSeqLength = std::vector< int >(1,300);
		ssIndInsModel = "PL";
		static double dIndLambdaDefault[2] = {0.0, 0.0};
		memcpy(dIndLambda, dIndLambdaDefault, sizeof(dIndLambda));
		dIndInsLambda = 0.0;
		bOutSeqGapSingleChar = false;
		dSubIota = std::vector< double >();
		dSubAlpha = std::vector< double >();
		dSeqReps = 1;
		ssOutBlockTail = "";
		ssOutFile = "-";
		dSubScale = std::vector< double >();
		nSubCodonWidth = 1;
		dSubGamma = std::vector< double >();
		ssIndModel[0] = "PL";
		ssIndModel[1] = "PL";
		dSeqRates = std::vector< std::vector< double > >();
		ssOutBlockAfter = "";
		static double dSubParamsDefault[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
		memcpy(dSubParams, dSubParamsDefault, sizeof(dSubParams));
		ssSubModel = "JC";
		ssSeqSequence = std::vector< std::string >(1,"");
		ssOutBlockBefore = "";

	}
	bool LoadStringOrFile(std::string &ss)
	{
		if(ss.empty())
			return false;
		std::ifstream iFile(ss.c_str());
		if(!iFile.is_open())
			return false;
		std::getline(iFile, ss, '\0');
		return true;
	}
	
};

} //namespace Dawg
#endif


