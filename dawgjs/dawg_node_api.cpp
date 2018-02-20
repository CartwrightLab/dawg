#include <string>
#include <random>
#include <sstream>

#include <dawg/subst.h>
#include <dawg/mutt.h>
#include <dawg/log.h>
#include <dawg/matic.h>

#include <boost/algorithm/string.hpp> // split

#include <emscripten/bind.h>

int getRandomInt(int low, int high) {
	std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist (low, high);
	return dist(mt);
}

float getRandomFloat(float low, float high) {
	std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist (low, high);
	return dist(mt);
}

class DawgWrapper {
public:
	DawgWrapper(unsigned int seed)
	: mRandomNumberGenerator() {
		mRandomNumberGenerator.seed(seed);
	}

	unsigned int getRandom(unsigned int a, unsigned int b) {
		using namespace std;
		// cout << "Hello rand(" << a << ", " << b << ")\n";
		auto n = mRandomNumberGenerator.rand_uint();
		// cout << "n: " << n << "\n";
		return n % b + a;
	}

	std::string getDawgError(const std::string &e) {
		using namespace std;
		stringstream ss;
		ss << DAWG_ERROR(e);
		return ss.str();
	}

	std::string getRnaSequence(const unsigned int root_length = 10000,
		const std::string &modelName = "jc") {

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

		// Model argument represents a section
		dawg::ma modelArgument;
		// initialize all variables directly unless they might be a list
		modelArgument.name = "argument1";
		// dawg::ma doesn't keep track of inheritance
		// modelArgument.inherits_from = inherits_from;

		modelArgument.subst_model = modelName;
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

		modelArgument.root_length = root_length;
		// modelArgument.root_seq = root_seq;
		// modelArgument.root_rates = splitIntoVectorDouble(root_rates);
		// modelArgument.root_code = root_code;
		// modelArgument.root_segment = root_segment;
		// modelArgument.root_gapoverlap = root_gapoverlap;

		modelArgument.output_rna = true;
		// modelArgument.output_keepempty = output_keepempty;
		modelArgument.output_markins = true;
		// modelArgument.output_lowercase = output_lowercase;

		std::vector<dawg::ma> modelArguments;
		modelArguments.emplace_back(modelArgument);

		Argument args (444, 10, "", "fasta:-");
		// for(string &ss : args.trick) {
		// WARNING: parse function has unresolve symbosl from EMCC compiler
		// ret &= input.parse(args.trick.begin(), args.trick.end());
		// }

		unsigned int num_reps = args.reps;

		dawg::matic kimura;

		kimura.seed(args.seed);

		if(!kimura.configure(modelArguments.begin(), modelArguments.end())) {
			DAWG_ERROR("bad configuration");
		}

		dawg::alignment aln;
		kimura.pre_walk(aln);
		for(unsigned int i=0;i<num_reps;++i) {
			kimura.walk(aln);
			mAlignments.emplace_back(aln);
			// write_aln(aln);
		}
		return getAlignments();
	}
private:
	dawg::mutt mRandomNumberGenerator;
	std::vector<dawg::alignment> mAlignments;

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

	std::string getAlignments() const {
		using namespace std;
		string braids;
		for (const auto &aln : mAlignments) {
			for (const auto &s : aln) {
				braids.append(s.label + ":" + s.seq + ";");
			}
		}
		// remove trailing semicolon
		braids.resize(braids.size() - 1);
		braids.shrink_to_fit();
		return braids;
	}
};

EMSCRIPTEN_BINDINGS(dawg_module) {
	using namespace emscripten;

	function("getRandomInt", &getRandomInt);
	function("getRandomFloat", &getRandomFloat);

	class_<DawgWrapper>("DawgWrapper")
	.constructor<unsigned int>()
	.function("getRandom", &DawgWrapper::getRandom)
	.function("getDawgError", &DawgWrapper::getDawgError)
	.function("getRnaSequence", &DawgWrapper::getRnaSequence)
	// .class_function("getStringFromInstance", &MyClass::getStringFromInstance)
	;
}
