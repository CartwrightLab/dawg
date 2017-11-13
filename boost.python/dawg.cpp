#include <vector>
#include <iostream>
#include <string>

#include <dawg/trick.h>
#include <dawg/trick_parse.h>
#include <dawg/global.h>
#include <dawg/output.h>
#include <dawg/ma.h>
#include <dawg/matic.h>

class DawgWalker {
public:
    void parse() {
        // std::vector<unsigned char> modifiedInstructionChars;
    	// input.parse(modifiedInstructionChars.begin(), modifiedInstructionChars.end());

    	// // process aliases
    	// input.read_aliases();
        //
    	// glopts.read_section(input.data.front());
        //
    	// // Since no output file has been specified,
    	// // the output will go to std::cout
    	// if (!write_aln.open(/*glopts.output_file.c_str()*/ "aln:-",
    	// 	glopts.sim_reps - 1,
    	// 	glopts.output_split, // false
    	// 	glopts.output_append, // false
    	// 	glopts.output_label)) // false
    	// {
    	// 	DAWG_ERROR("bad configuration");
    	// 	return;
    	// }
    	// write_aln.set_blocks(glopts.output_block_head.c_str(),
    	// 	glopts.output_block_between.c_str(),
    	// 	glopts.output_block_tail.c_str(),
    	// 	glopts.output_block_before.c_str(),
    	// 	glopts.output_block_after.c_str()
    	// );
    }

    void configure() {
        // std::vector<dawg::ma> configs;
    	// if (!dawg::ma::from_trick(input, configs)) {
    	// 	// DAWG_ERROR("bad configuration");
    	// 	return;
    	// }
        //
        // if (!kimura.configure(configs.begin(), configs.end())) {
        //     DAWG_ERROR("bad configuration");
        //     return;
        // }
    }

    void align() {
    	// create sets of aligned sequences;
    	// dawg::alignment aln;
    	// kimura.pre_walk(aln);
    	// for (unsigned int i = 0; i<glopts.sim_reps; ++i) {
    	// 	kimura.walk(aln);
    	// 	alignments.insert(alignments.end(), aln);
    	// 	//write_aln(aln); // this would print the aln data out to std::cout or a file
    	// }
    }

    void bark(const std::string& str) {
        std::cout << str << std::endl;
    }

private:
    // dawg::trick input;
    // dawg::matic kimura;
    // dawg::global_options glopts;
    // dawg::output write_aln;
    // std::vector<dawg::alignment> alignments;
};

#include <boost/python.hpp>

BOOST_PYTHON_MODULE(dawg_python)
{
    using namespace boost::python;
    class_<DawgWalker>("DawgWalker")
        .def("parse", &DawgWalker::parse)
        .def("configure", &DawgWalker::configure)
        .def("align", &DawgWalker::align)
        .def("bark", &DawgWalker::bark);
}
