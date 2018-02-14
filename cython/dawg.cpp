#include "dawg.hpp"

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <cstdlib> // atof

#include <boost/algorithm/string.hpp>

#include <dawg/ma.h>
#include <dawg/matic.h>
#include <dawg/trick_parse.h>
#include <dawg/global.h>
#include <dawg/log.h>

using namespace dawg;

// Example constructor to make sure things work
Dawg::Dawg() {
    info("no param constructor", __FILE__, __LINE__);
}

///////////////////////////////////////////////////////////
/// \brief constructor for using rng
///////////////////////////////////////////////////////////
Dawg::Dawg(const unsigned int seed)
: mSeed(seed)
, mRng() {
    mRng.seed(seed);
}

///////////////////////////////////////////////////////////
/// \param input can be optional
/// For instance, "dawg --seed=111 basic-dna.dawg -o fasta:-"
/// in the form of { in, out, r, s }
/// Test input to see if we need to use Dawg's Trick parser
///////////////////////////////////////////////////////////
Dawg::Dawg(const std::string& in,
    const std::string& out,
    const unsigned int seed,
    const unsigned int reps)
: mInFile(in)
, mOutFile(out)
, mSeed(seed)
, mRepetitions(reps)
, mTrickster()
, mKimura()
, mModelArguments()
, mWriter()
 {
	if (!in.empty())
		parseInput();

	// Create the object that will do all the simulation
	// work for us.  Configure its sections.
	// dawg::matic kimura;
    // if a seed was specified, use it
	if(mSeed != 0) {
		mKimura.seed(this->mSeed);
	}
}

void Dawg::configureMatic() {
	if (!mKimura.configure(mModelArguments.begin(), mModelArguments.end())) {
		DAWG_ERROR("bad configuration");
	}
}

void Dawg::parseInput() {
	// Parse the input file
	bool ret = true;
	auto pos = mInFile.rfind(".dawg");
	if (pos != std::string::npos)
		ret &= dawg::trick::parse_file(mTrickster, mInFile.c_str());
	else
		ret &= mTrickster.parse(mInFile.begin(), mInFile.end());

	if(!ret)
		DAWG_ERROR("Failure to parse DAWG file\n");

	// process aliases
	mTrickster.read_aliases();

	if (!dawg::ma::from_trick(mTrickster, mModelArguments)) {
		DAWG_ERROR("bad configuration: ");
	}
}

///////////////////////////////////////////////////////////
/// \brief addModelArgument
/// The traditional inheritance model in DAWG copies params
/// "downward" with respect to a trick file unless a section
///	specifies "[[secC = secA]]" the direction of inheritance
/// default params are set on Python side
///////////////////////////////////////////////////////////
void dawg::Dawg::addModelArgument(const std::string &name,
	const std::string &inherits_from,
	const std::string &subst_model,
    const std::string &subst_params,
    const std::string &subst_freqs,
    const std::string &subst_rate_model,
    const std::string &subst_rate_params,

    const std::string &indel_model_ins,
    const std::string &indel_params_ins,
    const std::string &indel_rate_ins,
    const unsigned int indel_max_ins,
    const std::string &indel_model_del,
    const std::string &indel_params_del,
    const std::string &indel_rate_del,
    const unsigned int indel_max_del,

	const std::string &tree,
    const std::string &tree_model,
    const std::string &tree_params,
    const double tree_scale,

    const unsigned int root_length,
    const std::string &root_seq,
    const std::string &root_rates,
    const unsigned int root_code,
    const unsigned int root_segment,
    const bool root_gapoverlap,

	const bool output_rna,
	const bool output_lowercase,
	const bool output_keepempty,
	const bool output_markins) {

	using namespace std;

	dawg::ma modelArgument;

	// initialize all variables directly unless they might be a list
	modelArgument.name = name;
	// dawg::ma doesn't keep track of inheritance
	// modelArgument.inherits_from = inherits_from;

	modelArgument.subst_model = subst_model;
	modelArgument.subst_params = splitIntoVectorDouble(subst_params);
	modelArgument.subst_freqs = splitIntoVectorDouble(subst_freqs);
	modelArgument.subst_rate_model = subst_rate_model;
	modelArgument.subst_rate_params = splitIntoVectorDouble(subst_rate_params);

	modelArgument.indel_model_ins = splitIntoVectorString(indel_model_ins);
	modelArgument.indel_params_ins = splitIntoVectorDouble(indel_params_ins);
	modelArgument.indel_rate_ins = splitIntoVectorDouble(indel_rate_ins);
	modelArgument.indel_max_ins = indel_max_ins;
	modelArgument.indel_model_del = splitIntoVectorString(indel_model_del);
	modelArgument.indel_params_del = splitIntoVectorDouble(indel_params_del);
	modelArgument.indel_rate_del = splitIntoVectorDouble(indel_rate_del);
	modelArgument.indel_max_del = indel_max_del;

	modelArgument.tree_tree = tree;
	modelArgument.tree_params = splitIntoVectorDouble(tree_params);
	modelArgument.tree_model = tree_model;
	modelArgument.tree_scale = tree_scale;

	modelArgument.root_length = root_length;
	modelArgument.root_seq = root_seq;
	// modelArgument.root_rates = splitIntoVectorDouble(root_rates);
	modelArgument.root_code = root_code;
	modelArgument.root_segment = root_segment;
	modelArgument.root_gapoverlap = root_gapoverlap;

	modelArgument.output_rna = output_rna;
	modelArgument.output_keepempty = output_keepempty;
	modelArgument.output_markins = output_markins;
	modelArgument.output_lowercase = output_lowercase;

	mModelArguments.emplace_back(modelArgument);
}

///////////////////////////////////////////////////////////
/// \brief Create the alignments
///	tell the dawg to walk, but then just walk it
///////////////////////////////////////////////////////////
void dawg::Dawg::walk()
{
    using namespace std;

	// prepare sets of aligned sequences;
	dawg::alignment aln;
	mKimura.pre_walk(aln);
	for (auto i = 0; i < mRepetitions; ++i) {
		mKimura.walk(aln);
		mAlignments.emplace_back(aln);
    }

} // run

///////////////////////////////////////////////////////////
/// \brief this would print the aln data out to std::cout
/// or a file depending on how output is defined
///////////////////////////////////////////////////////////
void dawg::Dawg::write() {
	using namespace std;

	dawg::global_options glopts;
	glopts.read_section(mTrickster.data.front());

	if (!mWriter.open(/*glopts.output_file.c_str()*/ mOutFile.c_str(),
		mRepetitions - 1,
		false, // false
		false, // false
		false)) // false
	{
		DAWG_ERROR("bad configuration");
		return;
	}

	for (auto a : mAlignments) {
		mWriter(a);
	}
}

///////////////////////////////////////////////////////////
/// \brief This is hacky and imperfect
/// We assume that the 'align' method in Dawg's evolver methods
/// are commented out. And then we just mash the alignment vector
///	together and return it as one giant string
///////////////////////////////////////////////////////////
std::string Dawg::getEvolvedSequences() const {
	using namespace std;
	string temp; // the string to append to
	for (const auto &aln : mAlignments) {
		for (const auto &s : aln) {
			temp += s.label + s.seq + ":";
		}
	}
	return temp;
}

unsigned int
dawg::Dawg::rand(unsigned int a, unsigned int b) {
    using namespace std;
    // cout << "Hello rand(" << a << ", " << b << ")\n";
    auto n = mRng.rand_uint();
    // cout << "n: " << n << "\n";
    return n % b + a;
}

///////////////////////////////////////////////////////////
///	\brief bark :: Print out segments (model args) to stdout
///////////////////////////////////////////////////////////
void Dawg::bark() const {
	using namespace std;

	for (auto arg : mModelArguments) {
		cout << arg << endl;
	}
}

std::vector<std::string> Dawg::splitIntoVectorString(const std::string &s) const {
	using namespace std;

	vector<string> string_split;
	boost::algorithm::split(string_split, s, boost::is_any_of(","));

	return string_split;
}

std::vector<double> Dawg::splitIntoVectorDouble(const std::string &s) const {
	using namespace std;

	vector<string> string_split;
	boost::algorithm::split(string_split, s, boost::is_any_of(","));

	vector<double> string_to_double;
	for (auto params : string_split) {
		string_to_double.emplace_back(atof(params.c_str())); // consider strtod
	}
	return string_to_double;
}

void Dawg::printAlignment(const dawg::alignment &aln) const {
	using namespace std;
	// cout << "mAlignments.size(): " << mAlignments.size() << endl;
	cout << "max_label_width: " << aln.max_label_width <<
		"\nseq_type: " << aln.seq_type << "\n";
	for (const auto &v : aln) {
		cout << "label: " << v.label.c_str() << "\nseq: " << v.seq.c_str() << endl;
	}
}

// Log an info message
template <typename Line, typename File>
void Dawg::info(const std::string &msg, Line l, File f) const {
	using namespace std;
	cout << msg << ", file: " << f << ", line: " << l << endl;
}

// Log an error message
template <typename Line, typename File>
void Dawg::error(const std::string &msg, Line l, File f) const {
	using namespace std;
	cerr << msg << ", file: " << f << ", line: " << l << endl;
}
