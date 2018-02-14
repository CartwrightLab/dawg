#ifndef DAWG_HPP
#define DAWG_HPP

#include <string>
#include <vector>

#include <dawg/mutt.h>
#include <dawg/trick.h>
#include <dawg/residue.h>
#include <dawg/matic.h>
#include <dawg/ma.h>
#include <dawg/output.h>

namespace dawg {

class Dawg
{
public:
    explicit Dawg();
    explicit Dawg(const unsigned int seed);
    explicit Dawg(const std::string& in,
        const std::string& out,
        const unsigned int seed,
        const unsigned int reps);

    void addModelArgument(const std::string &name,
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
        const bool output_markins);
    void configureMatic();
    void walk();
    void write();
    std::string getEvolvedSequences() const;
    unsigned int rand(unsigned int a, unsigned int b);
    void bark() const;
private:
    std::string mInFile;
    std::string mOutFile;
    unsigned int mSeed;
    std::size_t mRepetitions;
    std::vector<dawg::alignment> mAlignments;
    mutt mRng;
    trick mTrickster;
    dawg::matic mKimura;
    std::vector<dawg::ma> mModelArguments;
    dawg::output mWriter;

    void parseInput();

    std::vector<std::string> splitIntoVectorString(const std::string &s) const;
    std::vector<double> splitIntoVectorDouble(const std::string &s) const;

    void printAlignment(const dawg::alignment &aln) const;

    template <typename Line, typename File>
    void info(const std::string &msg, Line l, File f) const;
    template <typename Line, typename File>
    void error(const std::string &msg, Line l, File f) const;
}; // class Dawg

} // namespace dawg

#endif
