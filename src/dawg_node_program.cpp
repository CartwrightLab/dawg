/*  Dawg - DNA Assembly with Gaps - Simulating Sequence Evolution
    Copyright (c) 2004-2012  Reed A. Cartwright, PhD

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
#include "../build/release/version.h" // hack for this branch
// #include "dawg.h"

#include <boost/preprocessor.hpp>
// #include <boost/config.hpp>

#include <vector>
#include <exception>

#include <dawg/matic.h>
#include <dawg/ma.h>
#include <dawg/trick.h>
#include <dawg/global.h>
#include <dawg/output.h>
#include <dawg/log.h>

// #include "dawg_app.h"

using namespace std;
using namespace boost;
using namespace dawg;

#define VERSION_MSG NEW_PACKAGE_STRING "\n" \
	"    Copyright (C) 2004-2013  Reed A. Cartwright, PhD <cartwright@asu.edu>\n"

std::vector<std::string> splitIntoVectorString(const std::string &s);
std::vector<double> splitIntoVectorDouble(const std::string &s);
/// \brief Run the application
int run();

int main(int argc, char *argv[])
{
	int ret = EXIT_FAILURE;
	try {
		// dawg_app app(argc, argv);
		// ret = app.run();
		run();
	} catch(std::exception &e) {
		DAWG_ERROR(e.what());
	}
	return ret;
}

std::vector<std::string> splitIntoVectorString(const std::string &s) {
	using namespace std;

	vector<string> string_split;
	boost::algorithm::split(string_split, s, boost::is_any_of(","));

	return string_split;
}

std::vector<double> splitIntoVectorDouble(const std::string &s) {
	using namespace std;

	vector<string> string_split;
	boost::algorithm::split(string_split, s, boost::is_any_of(","));

	vector<double> string_to_double;
	for (auto params : string_split) {
		string_to_double.emplace_back(atof(params.c_str())); // consider strtod
	}
	return string_to_double;
}

int run() {

	class Argument {
	public:
		unsigned int seed, reps;
		std::string trick, output;

		explicit Argument(
			unsigned int seed,
			unsigned int reps,
			const std::string &trick,
			const std::string &output)
		: seed(seed)
		, reps(reps)
		, trick(trick)
		, output(output) {

		}
	};

	std::string trickString = " \
		[[-]] \
		Subst.Model = jc \
		Root.Segment = 1 \
		Root.Length = 60 \
		Tree.Tree = ((A:0.02,B:0.02):0.2,(C:0.02):0.2); \
		[[-]] \
		Root.Code = 1 \
		Root.Segment = 0 \
		Root.Length = 60 \
		Tree.Tree = ((A:0.02):0.2,(B:0.02,C:0.02):0.2); \
		[[-]] \
		Root.Segment = 2";

		// Model argument represents a section
		dawg::ma modelArgument;
		// initialize all variables directly unless they might be a list
		modelArgument.name = "argument1";
		// dawg::ma doesn't keep track of inheritance
		// modelArgument.inherits_from = inherits_from;

		modelArgument.subst_model = "jc";
		// modelArgument.subst_params = splitIntoVectorDouble(subst_params);
		// modelArgument.subst_freqs = splitIntoVectorDouble(subst_freqs);
		// modelArgument.subst_rate_model = subst_rate_model;
		// modelArgument.subst_rate_params = splitIntoVectorDouble(subst_rate_params);

		// modelArgument.indel_model_ins = splitIntoVectorString(indel_model_ins);
		// modelArgument.indel_params_ins = splitIntoVectorDouble(indel_params_ins);
		// modelArgument.indel_rate_ins = splitIntoVectorDouble(indel_rate_ins);
		// modelArgument.indel_max_ins = indel_max_ins;
		// modelArgument.indel_model_del = splitIntoVectorString(indel_model_del);
		// modelArgument.indel_params_del = splitIntoVectorDouble(indel_params_del);
		// modelArgument.indel_rate_del = splitIntoVectorDouble(indel_rate_del);
		// modelArgument.indel_max_del = indel_max_del;

		modelArgument.tree_tree = "((A:0.02,B:0.02):0.2,(C:0.02):0.2);";
		// modelArgument.tree_params = splitIntoVectorDouble(tree_params);
		// modelArgument.tree_model = tree_model;
		// modelArgument.tree_scale = tree_scale;

		modelArgument.root_length = 60;
		// modelArgument.root_seq = root_seq;
		// modelArgument.root_rates = splitIntoVectorDouble(root_rates);
		// modelArgument.root_code = root_code;
		// modelArgument.root_segment = root_segment;
		// modelArgument.root_gapoverlap = root_gapoverlap;

		// modelArgument.output_rna = output_rna;
		// modelArgument.output_keepempty = output_keepempty;
		// modelArgument.output_markins = output_markins;
		// modelArgument.output_lowercase = output_lowercase;

	std::vector<dawg::ma> modelArguments;
	modelArguments.emplace_back(modelArgument);

	Argument args (444, 10, trickString, "fasta:-");


	//if(arg.quiet)
	//	cerr.clear(ios::failbit);
	trick input;

	bool ret = true;
	// for(string &ss : args.trick) {
	// WARNING: parse function has unresolve symbosl from EMCC compiler
	// ret &= input.parse(args.trick.begin(), args.trick.end());
	// }

	if(!ret)
		return EXIT_FAILURE;
	// process aliases
	input.read_aliases();

	global_options glopts;
	glopts.read_section(input.data.front());

	unsigned int num_reps = args.reps; // > 0) ? args.reps : glopts.sim_reps;

	dawg::output write_aln;
	const char *file_name = args.output.c_str(); //arg.output.empty() ? glopts.output_file.c_str() : arg.output.c_str();
	//bool split  = (!vm["split"].defaulted()) ? arg.split : glopts.output_split;
	//bool append = (!vm["append"].defaulted()) ? arg.append : glopts.output_append;

	bool split  = false ; //arg.split || (indeterminate(arg.split) && glopts.output_split);
	bool append = false; //arg.append || (indeterminate(arg.append) && glopts.output_append);
	bool label  = false ; //arg.label || (indeterminate(arg.label) && glopts.output_label);

	if(!write_aln.open(file_name, num_reps-1, split, append, label)) {
		DAWG_ERROR("bad configuration");
		return EXIT_FAILURE;
	}

	vector<dawg::ma> configs;
	if(!dawg::ma::from_trick(input, configs)) {
		DAWG_ERROR("bad configuration");
		return EXIT_FAILURE;
	}

	// Create the object that will do all the simulation
	// work for us.  Configure its sections.
	dawg::matic kimura;
	// if a seed was specified, use it
	// if(arg.seed != 0) {
	kimura.seed(args.seed);
	// } else if(!glopts.sim_seed.empty()) {
		// kimura.seed(glopts.sim_seed.begin(), glopts.sim_seed.end());
	// }
	if(!kimura.configure(modelArguments.begin(), modelArguments.end())) {
		DAWG_ERROR("bad configuration");
		return EXIT_FAILURE;
	}
	// create sets of aligned sequences;
	dawg::alignment aln;
	kimura.pre_walk(aln);
	for(unsigned int i=0;i<num_reps;++i) {
		kimura.walk(aln);
		write_aln(aln);
	}
	return EXIT_SUCCESS;
} // run
