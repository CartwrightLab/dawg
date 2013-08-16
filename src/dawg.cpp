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
#include "../version.h"
#include "dawg.h"

#include <dawg/matic.h>

#include <boost/preprocessor.hpp>
#include <boost/config.hpp>

#include <exception>

#include <dawg/ma.h>
#include <dawg/trick.h>
#include <dawg/utils/foreach.h>
#include <dawg/global.h>
#include <dawg/output.h>

#include "dawg_app.h"

using namespace std;
using namespace boost;
using namespace dawg;

#define VERSION_MSG NEW_PACKAGE_STRING "\n" \
	"    Copyright (C) 2004-2013  Reed A. Cartwright, PhD <cartwright@asu.edu>\n"

int main(int argc, char *argv[])
{
	int ret = EXIT_FAILURE;
	try {
		dawg_app app(argc, argv);
		ret = app.run();
	} catch(std::exception &e) {
		CERROR(e.what());
	}
	return ret;
}

dawg_app::dawg_app(int argc, char* argv[]) : desc("Allowed Options") {
	runname = argv[0];
	try {
		desc.add_options()
			#define XM(lname, sname, desc, type, def) ( \
				XS(lname) IFD(sname, "," BOOST_PP_STRINGIZE sname), \
				po::value< type >(&arg.XV(lname))->default_value(def), \
				desc )				
			#include "dawgarg.xmh"
			#undef XM
			;
		indesc.add_options()("input", po::value< vector<string> >(&arg.input), "input files");
		indesc.add(desc);
		pdesc.add("input", -1);
		po::store(po::command_line_parser(argc, argv).options(indesc).positional(pdesc).run(), vm);
		po::notify(vm);
		if(!arg.arg_file.empty()) {
			if(arg.arg_file == "-") {
				po::store(po::parse_config_file(cin, desc), vm);	
			} else {
				std::ifstream ifs(arg.arg_file.c_str());
				if(!ifs.is_open()) {
					string sse = "unable to open argument file ";
					sse += arg.arg_file;
					throw std::runtime_error(sse);
				}
				po::store(po::parse_config_file(ifs, desc), vm);
			}
			po::notify(vm);
		}
	} catch (std::exception &e) {
		CERROR(e.what());
		throw std::runtime_error("unable to process command line");
	}
}

int dawg_app::run() {
	//std::string _temp(" I love %r/%R/%%/%.  Do you?"), _out;
	//_out = boost::algorithm::replace_all_regex_copy(_temp, boost::regex("%(r)|%(R)|%(%)"),
	//	std::string("?1x:?2y:z"), match_default | format_all);
	//cout << _out << endl << endl;
	

	if(arg.version)	{
		cerr << endl << VERSION_MSG << endl << endl;
		return EXIT_SUCCESS;
	}
	if(arg.help_trick) {
		cerr << endl << VERSION_MSG << endl << endl;
		ma::help(cerr);
		return EXIT_SUCCESS;
	}
	if(arg.help || arg.input.empty()) {
		cerr << endl << VERSION_MSG << endl << endl;
		cerr << "Usage:\n  "
		     << runname << " [options] trick-1.dawg trick-2.dawg ..."
			 << endl << endl;
		cerr << desc << endl;
		return EXIT_SUCCESS;
	}
	
	//if(arg.quiet)
	//	cerr.clear(ios::failbit);
	trick input;
		
	bool ret = true;
	foreach(string &ss, arg.input) {
		ret &= trick::parse_file(input, ss.c_str());
	}
	
	if(!ret)
		return EXIT_FAILURE;	
	// process aliases
	input.read_aliases();

	global_options glopts;
	glopts.read_section(input.data.front());

	unsigned int num_reps = (arg.reps > 0) ? arg.reps : glopts.sim_reps;

	dawg::output write_aln;
	const char *file_name = arg.output.empty() ? glopts.output_file.c_str()
		                                       : arg.output.c_str();
	//bool split  = (!vm["split"].defaulted()) ? arg.split : glopts.output_split;
	//bool append = (!vm["append"].defaulted()) ? arg.append : glopts.output_append;

	bool split  = arg.split || (indeterminate(arg.split) && glopts.output_split);
	bool append = arg.append || (indeterminate(arg.append) && glopts.output_append);
	bool label  = arg.label || (indeterminate(arg.label) && glopts.output_label);

	if(!write_aln.open(file_name, num_reps-1, split, append, label)) {
		DAWG_ERROR("bad configuration");
		return EXIT_FAILURE;
	}
	write_aln.set_blocks(glopts.output_block_head.c_str(),
		glopts.output_block_between.c_str(),
		glopts.output_block_tail.c_str(),
		glopts.output_block_before.c_str(),
		glopts.output_block_after.c_str()
	);

	vector<dawg::ma> configs;
	if(!dawg::ma::from_trick(input, configs)) {
		DAWG_ERROR("bad configuration");
		return EXIT_FAILURE;
	}

	// Create the object that will do all the simulation
	// work for us.  Configure its sections.
	dawg::matic kimura;
	// if a seed was specified, use it
	if(arg.seed != 0) {
		kimura.seed(arg.seed);
	} else if(!glopts.sim_seed.empty()) {
		kimura.seed(glopts.sim_seed.begin(), glopts.sim_seed.end());
	}
	if(!kimura.configure(configs.begin(), configs.end())) {
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
}
