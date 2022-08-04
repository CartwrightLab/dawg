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
#include "dawg.h"

#include <CLI11.hpp>
#include <boost/preprocessor.hpp>
#include <exception>

#include "../version.h"
#include "dawg/global.h"
#include "dawg/ma.h"
#include "dawg/matic.h"
#include "dawg/output.h"
#include "dawg/trick.h"
#include "dawg_app.h"

#define VERSION_MSG                                         \
    DAWG_PACKAGE_STRING                                     \
    "\n"                                                    \
    "    Copyright (C) 2004-2013  Reed A. Cartwright, PhD " \
    "<cartwright@asu.edu>\n"

int main(int argc, char *argv[]) {
    int ret = EXIT_FAILURE;
    try {
        dawg_app app(argc, argv);
        ret = app.run();
    } catch(std::exception &e) {
        CERROR(e.what());
    }
    return ret;
}

dawg_app::dawg_app(int argc, char *argv[]) : runname(argv[0]) {
    // set_cli_options
    this->cli_app.add_option("input", arg.input, "input files");
#define XM(lname, sname, desc, type, def)                            \
    this->cli_app.add_option(                                        \
        IFD(sname, "-" BOOST_PP_STRINGIZE sname ",") "--" XS(lname), \
        arg.XV(lname), desc, def);
#define XF(lname, sname, desc, type, def)                            \
    this->cli_app.add_flag(                                          \
        IFD(sname, "-" BOOST_PP_STRINGIZE sname ",") "--" XS(lname), \
        arg.XV(lname), desc);
#include "dawgarg.xmh"
#undef XM
#undef XF

    try {
        this->cli_app.parse(argc, argv);
    } catch(const CLI::CallForHelp &e) {
        exit(this->cli_app.exit(e));
    }
}

int dawg_app::run() {
    // std::string _temp(" I love %r/%R/%%/%.  Do you?"), _out;
    //_out = boost::algorithm::replace_all_regex_copy(_temp,
    // boost::regex("%(r)|%(R)|%(%)"), 	std::string("?1x:?2y:z"), match_default
    // | format_all); cout << _out << endl << endl;

    if(arg.version) {
        std::cerr << std::endl << VERSION_MSG << std::endl << std::endl;
        return EXIT_SUCCESS;
    }
    if(arg.help_trick) {
        std::cerr << std::endl << VERSION_MSG << std::endl << std::endl;
        dawg::ma::help(std::cerr);
        return EXIT_SUCCESS;
    }
    if(arg.input.empty()) {
        std::cerr << std::endl << VERSION_MSG << std::endl << std::endl;
        std::cerr << std::endl << this->cli_app.help() << std::endl;
        return EXIT_SUCCESS;
    }

    // if(arg.quiet)
    //	cerr.clear(ios::failbit);
    dawg::trick input;

    bool ret = true;
    for(std::string &ss : arg.input) {
        ret &= dawg::trick::parse_file(input, ss.c_str());
    }

    if(!ret) return EXIT_FAILURE;
    // process aliases
    input.read_aliases();

    dawg::global_options glopts;
    glopts.read_section(input.data.front());

    unsigned int num_reps = (arg.reps > 0) ? arg.reps : glopts.sim_reps;

    dawg::output write_aln;
    const char *file_name =
        arg.output.empty() ? glopts.output_file.c_str() : arg.output.c_str();

    if(!write_aln.open(file_name, num_reps - 1, arg.split, arg.append,
                       arg.label)) {
        DAWG_ERROR("bad configuration");
        return EXIT_FAILURE;
    }
    write_aln.set_blocks(
        glopts.output_block_head.c_str(), glopts.output_block_between.c_str(),
        glopts.output_block_tail.c_str(), glopts.output_block_before.c_str(),
        glopts.output_block_after.c_str());

    std::vector<dawg::ma> configs;
    if(!dawg::ma::from_trick(input, configs)) {
        DAWG_ERROR("bad configuration");
        return EXIT_FAILURE;
    }

    // Create the object that will do all the simulation
    // work for us.  Configure its sections.
    dawg::matic simulation;
    // if a seed was specified, use it
    if(arg.seed != 0) {
        simulation.seed(arg.seed);
    } else if(!glopts.sim_seed.empty()) {
        simulation.seed(glopts.sim_seed.begin(), glopts.sim_seed.end());
    } else {
        simulation.auto_seed_seq();
    }

    if(!simulation.configure(configs.begin(), configs.end())) {
        DAWG_ERROR("bad configuration");
        return EXIT_FAILURE;
    }
    // create sets of aligned sequences;
    dawg::alignment aln;
    simulation.pre_walk(aln);
    for(unsigned int i = 0; i < num_reps; ++i) {
        simulation.walk(aln);
        write_aln(aln);
    }
    return EXIT_SUCCESS;
}
