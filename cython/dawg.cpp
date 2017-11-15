#include "dawg.hpp"

#include "../src/dawg_app.h"
#include <dawg/ma.h>
#include <dawg/trick.h>
#include <dawg/global.h>
#include <dawg/output.h>

#include <iostream>

dawg::Dawg::Dawg()
{
    using namespace std;
    cout << "Hello PyDawg: " << __FILE__ << ", " << __LINE__ << endl;
}

dawg::Dawg::Dawg(const std::string& s)
{
    using namespace std;
    cout << "Hello PyDawg: " << __FILE__ << ", " << __LINE__ << endl;
    cout << "s: " << s << "\n";
}

void dawg::Dawg::run()
{
/*
#ifdef WITH_DAWG
vector<string> Simulation::tempDawgStrToParamSpecDawgStr(int numberOfSimulations, const vector<string> & templateInstructionString)
{
	vector<string> modifiedInstructionString;
	for (int i = 0; i<templateInstructionString.size(); i++)
	{
		string line = templateInstructionString[i];
		if (line.find("Params.Ins") != std::string::npos)
		{
			string newLine = "Params.Ins = " + to_string(static_cast<long double>(_indelDistributionShapeParameter)) + ", " + to_string(static_cast<long long>(_MaxdeletionLengh)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Params.Del") != std::string::npos)
		{
			string newLine = "Params.Del = " + to_string(static_cast<long double>(_indelDistributionShapeParameter)) + ", " + to_string(static_cast<long long>(_MaxdeletionLengh)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Rate.Ins") != std::string::npos)
		{
			string newLine = "Rate.Ins = " + to_string(static_cast<long double>(_indelRate)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Rate.Del") != std::string::npos)
		{
			string newLine = "Rate.Del = " + to_string(static_cast<long double>(_indelRate)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else if (line.find("Length") != std::string::npos)
		{
			string newLine = "Length = " + to_string(static_cast<long long>(_rootLength)) + "\n";
			modifiedInstructionString.push_back(newLine);
		}
		else
		{
			modifiedInstructionString.push_back(line);
		}
	}
	return modifiedInstructionString;
}

void Simulation::simulateDawgMSA(const int numberOfSimulations, const vector<string> & templateInstructionString)
{
	vector<string> modifiedInstructionString = tempDawgStrToParamSpecDawgStr(numberOfSimulations, templateInstructionString);
	vector<char> modifiedInstructionChars = convertVecStringToVecChar(modifiedInstructionString);
	generateDawgMSA_array(numberOfSimulations, modifiedInstructionChars);
}

void Simulation::generateDawgMSA_array(int numberOfSimulations, const vector<char> & modifiedInstructionChars)
{

	msaVec.clear();
	// The alignments created from DAWG
    std::vector<dawg::alignment> alignments;
	for (int i = 0; i < numberOfSimulations; ++i)
	{
		createDawgAlignments(alignments, modifiedInstructionChars);
		dawg::alignment single_dawg_alignment = alignments[0];
		vector<string> dawg_aligned_seqs;
		for (int i = 0; i < single_dawg_alignment.size(); i++)
		{
			//string name = single_dawg_alignment[i].label;
			//string seq = single_dawg_alignment[i].seq;
			//cout << name << "\t" << seq << endl;

			dawg_aligned_seqs.push_back(single_dawg_alignment[i].seq);
		}

		MSA dawg_sim_MSA(dawg_aligned_seqs);
		msaVec.push_back(dawg_sim_MSA);
		alignments.clear();
	}
}*/



// Create alignments

	dawg::trick input;

	input.parse(modifiedInstructionChars.begin(), modifiedInstructionChars.end());

	// process aliases
	input.read_aliases();

	dawg::global_options glopts;
	glopts.read_section(input.data.front());

	dawg::output write_aln;

	// Since no output file has been specified,
	// the output will go to std::cout
	if (!write_aln.open(/*glopts.output_file.c_str()*/ "aln:-",
		glopts.sim_reps - 1,
		glopts.output_split, // false
		glopts.output_append, // false
		glopts.output_label)) // false
	{
		DAWG_ERROR("bad configuration");
		return;
	}
	write_aln.set_blocks(glopts.output_block_head.c_str(),
		glopts.output_block_between.c_str(),
		glopts.output_block_tail.c_str(),
		glopts.output_block_before.c_str(),
		glopts.output_block_after.c_str()
	);

	std::vector<dawg::ma> configs;
	if (!dawg::ma::from_trick(input, configs)) {
		DAWG_ERROR("bad configuration");
		return;
	}

	// Create the object that will do all the simulation
	// work for us.  Configure its sections.
	dawg::matic kimura;

	if (!kimura.configure(configs.begin(), configs.end())) {
		DAWG_ERROR("bad configuration");
		return;
	}

	// create sets of aligned sequences;
	dawg::alignment aln;
	kimura.pre_walk(aln);
	for (unsigned int i = 0; i<glopts.sim_reps; ++i) {
		kimura.walk(aln);
		alignments.insert(alignments.end(), aln);
		//write_aln(aln); // this would print the aln data out to std::cout or a file
}

}
