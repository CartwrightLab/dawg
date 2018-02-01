#include "dawg.hpp"

#include <iostream>
#include <cstring>
#include <string>
#include <vector>

#include <dawg/ma.h>
#include <dawg/matic.h>
#include <dawg/trick_parse.h>
#include <dawg/global.h>
#include <dawg/output.h>

dawg::Dawg::Dawg()
{
    using namespace std;
    cout << "Hello PyDawg: " << __FILE__ << ", " << __LINE__ << endl;
}

dawg::Dawg::Dawg(const unsigned int s)
: seed(s)
, rng()
{
    using namespace std;
    rng.seed(s);
}

///////////////////////////////////////////////////////////
/// \param s is the input string
/// For instance, "dawg --seed=111 basic-dna.dawg -o fasta:-"
/// in the form of { in, out, r, s }
///////////////////////////////////////////////////////////
dawg::Dawg::Dawg(const std::string& in,
    const std::string& out,
    const unsigned int r,
    const unsigned int s)
: inFile(in)
, outFile(out)
, reps(r)
, seed(s)
, mAlignments(reps)
{

}

///////////////////////////////////////////////////////////
/// \brief Create the alignments
///
///////////////////////////////////////////////////////////
void dawg::Dawg::run()
{
    using namespace std;

    bool ret = true;
	// for(string &ss : inFile) {
		ret &= dawg::trick::parse_file(mInput, inFile.c_str());
	// }

	if(!ret)
		std::cerr << "Failure to parse DAWG file\n";

	// process aliases
	mInput.read_aliases();

	dawg::global_options glopts;
	glopts.read_section(mInput.data.front());

	dawg::output write_aln;

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
	if (!dawg::ma::from_trick(mInput, configs)) { // throwing error
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

	// prepare sets of aligned sequences;
	kimura.pre_walk(mAlignments.front());
	cout << "mAlignments.size(): " << mAlignments.size() << endl;
	for (auto i = 0; i < reps; ++i) {
		kimura.walk(mAlignments.at(i));
		// alignments.insert(alignments.end(), aln);
		// write_aln(mAlignments.at(i)); // this would print the aln data out to std::cout or a file
		printAlignmentInfo(mAlignments.at(i));
    }

} // run

void dawg::Dawg::bark() const {
    using namespace std;
    cout << "inFile: " << inFile << ", " <<
        "outFile: " << outFile << ", " <<
        "reps: " << reps << ", " <<
        "seed: " << seed << "\n";
} // bark

unsigned int
dawg::Dawg::rand(unsigned int a, unsigned int b) {
    using namespace std;
    // cout << "Hello rand(" << a << ", " << b << ")\n";
    auto n = rng.rand_uint();
    // cout << "n: " << n << "\n";
    return n % b + a;
}

void dawg::Dawg::trickStats() const {
    using namespace dawg;
    using namespace std;

    auto sections = mInput.data;
    for (auto sec : sections) {
        cout << "section: " << sec.name << "\n" <<
        "inherits: " << sec.inherits << "\n";
        for (auto node : sec.db) {
            cout << "node: " << node.first << ", values: ";
            for (auto value : node.second) {
                cout << value << ", ";
            }
            cout << "\n";
        }
    }
}

void dawg::Dawg::printAlignmentInfo(const dawg::alignment &aln) const {
	using namespace std;
	cout << "max_label_width: " << aln.max_label_width <<
		"\nseq_type: " << aln.seq_type << "\n";
	for (const auto &v : aln) {
		cout << "label: " << v.label.c_str() << "\nseq: " << v.seq.c_str() << endl;
	}
}