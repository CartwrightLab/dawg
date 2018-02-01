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

// Example constructor to make sure things work
dawg::Dawg::Dawg()
{
    using namespace std;
    cout << "Hello PyDawg: " << __FILE__ << ", " << __LINE__ << endl;
}

dawg::Dawg::Dawg(const unsigned int s)
: mSeed(s)
, mRng()
{
    using namespace std;
    mRng.seed(s);
}

dawg::Dawg::Dawg(const std::map<std::string, std::vector<std::string>>,
	const std::string& o,
	const unsigned int seed,
	const unsigned int reps)
{
    using namespace std;
    cout << "Hello PyDawg: " << __FILE__ << ", " << __LINE__ << endl;
}

///////////////////////////////////////////////////////////
/// \param s is the input string
/// For instance, "dawg --seed=111 basic-dna.dawg -o fasta:-"
/// in the form of { in, out, r, s }
///////////////////////////////////////////////////////////
dawg::Dawg::Dawg(const std::string& in,
    const std::string& out,
    const unsigned int seed,
    const unsigned int reps)
: mInFile(in)
, mOutFile(out)
, mSeed(seed)
, mRepetitions(reps)
, mTrickster()
{
	bool ret = true;

	auto pos = in.rfind(".dawg");

	if (pos != std::string::npos)
		ret &= dawg::trick::parse_file(mTrickster, mInFile.c_str());
	else
		ret &= mTrickster.parse(mInFile.begin(), mInFile.end());

	if(!ret)
		std::cerr << "Failure to parse DAWG file\n";
}

///////////////////////////////////////////////////////////
/// \brief Create the alignments
///
///////////////////////////////////////////////////////////
void dawg::Dawg::run()
{
    using namespace std;

	// process aliases
	mTrickster.read_aliases();

	std::vector<dawg::ma> configs;
	if (!dawg::ma::from_trick(mTrickster, configs)) {
		DAWG_ERROR("bad configuration");
		return;
	}

	// Create the object that will do all the simulation
	// work for us.  Configure its sections.
	dawg::matic kimura;
    // if a seed was specified, use it
	if(mSeed != 0) {
		kimura.seed(this->mSeed);
	}

	if (!kimura.configure(configs.begin(), configs.end())) {
		DAWG_ERROR("bad configuration");
		return;
	}

	// prepare sets of aligned sequences;
	dawg::alignment aln;
	kimura.pre_walk(aln);
	for (auto i = 0; i < mRepetitions; ++i) {
		kimura.walk(aln);
		// printAlignmentInfo(aln);
		// mAlignments.at(i) = std::move(aln);
		mAlignments.emplace_back(aln);
    }

} // run

///////////////////////////////////////////////////////////
/// \brief this would print the aln data out to std::cout or a file
///////////////////////////////////////////////////////////
void dawg::Dawg::printAlignments() {
	using namespace std;
	dawg::output write_aln;

	dawg::global_options glopts;
	glopts.read_section(mTrickster.data.front());

	if (!write_aln.open(/*glopts.output_file.c_str()*/ mOutFile.c_str(),
		mRepetitions - 1,
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

	for (auto a : mAlignments) {
#if defined(DAWG_DEBUG)
	printAlignmentInfo(a);
#endif // defined
	
		write_aln(a);
	}
}

void dawg::Dawg::bark() const {
    using namespace std;
    cout << "inFile: " << mInFile << ", " <<
        "outFile: " << mOutFile << ", " <<
        "reps: " << mRepetitions << ", " <<
        "seed: " << mSeed << "\n";
} // bark

unsigned int
dawg::Dawg::rand(unsigned int a, unsigned int b) {
    using namespace std;
    // cout << "Hello rand(" << a << ", " << b << ")\n";
    auto n = mRng.rand_uint();
    // cout << "n: " << n << "\n";
    return n % b + a;
}

void dawg::Dawg::trickStats() const {
    using namespace dawg;
    using namespace std;

    auto sections = mTrickster.data;
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
	// cout << "mAlignments.size(): " << mAlignments.size() << endl;
	cout << "max_label_width: " << aln.max_label_width <<
		"\nseq_type: " << aln.seq_type << "\n";
	for (const auto &v : aln) {
		cout << "label: " << v.label.c_str() << "\nseq: " << v.seq.c_str() << endl;
	}
}
