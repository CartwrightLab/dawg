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
    DawgWalker();
    DawgWalker(const std::string &in,
        const std::string &out,
        const unsigned int r,
        const unsigned int s)
    : inFile(in), outFile(out),
    reps(r), seed(s)
    {

    }

    void run() {
        using std::string;
        using dawg::trick;

    	dawg::trick input;

        bool ret = true;
    	// for(string &ss : inFile) {
    		ret &= trick::parse_file(input, inFile.c_str());
    	// }

    	if(!ret)
    		std::cerr << "Failure to parse DAWG file\n";

    	// process aliases
    	input.read_aliases();

    	dawg::global_options glopts;
    	glopts.read_section(input.data.front());

    	dawg::output write_aln;

    	// Since no output file has been specified,
    	// the output will go to std::cout
    	if (!write_aln.open(/*glopts.output_file.c_str()*/ outFile.c_str(),
    		reps - 1,
    		false, // false
    		false, // false
    		false)) // false
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
        // if a seed was specified, use it
    	if(seed != 0) {
    		kimura.seed(this->seed);
    	}

    	if (!kimura.configure(configs.begin(), configs.end())) {
    		DAWG_ERROR("bad configuration");
    		return;
    	}

    	// create sets of aligned sequences;
        std::vector<dawg::alignment> alignments;
    	dawg::alignment aln;
    	kimura.pre_walk(aln);
    	for (unsigned int i = 0; i< reps; ++i) {
    		kimura.walk(aln);
    		alignments.insert(alignments.end(), aln);
    		//write_aln(aln); // this would print the aln data out to std::cout or a file
        }
    } // run

    void bark(const std::string &str) {
        std::cout << str << std::endl;
    }

private:
    const std::string &inFile, &outFile;
    const unsigned int reps, seed;
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
        .def("run", &DawgWalker::run)
        .def("bark", &DawgWalker::bark);
}
