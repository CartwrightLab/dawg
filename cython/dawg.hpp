#ifndef DAWG_HPP
#define DAWG_HPP

#include <string>
#include <vector>
#include <map>

#include <dawg/mutt.h>
#include <dawg/trick.h>
#include <dawg/residue.h>
#include <dawg/matic.h>
#include <dawg/ma.h>

namespace dawg {

class Dawg
{
public:
    explicit Dawg();
    explicit Dawg(const unsigned int s);

    explicit Dawg(const std::string& output,
        const unsigned int seed,
        const unsigned int reps);

    explicit Dawg(const std::string& infile,
        const std::string& outfile,
        const unsigned int seed,
        const unsigned int reps);

    void addSegment(const std::string &name,
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

        const std::string &tree_model,
        const std::string &tree_params,
        const std::string &tree_tree,
        const double tree_scale,

        const unsigned int root_length,
        const std::string &root_seq,
        const std::string &root_rates,
        const unsigned int root_code,
        const unsigned int root_segment,
        const bool root_gapoverlap,

        const bool output_markins,
        const bool output_keepempty,
        const bool output_lowercase,
        const bool output_rna);

    void echoSegments() const;

    void run();
    void printAlignments();
    void bark() const;
    unsigned int rand(unsigned int a, unsigned int b);
    void trickStats() const;

    const char * getEvolvedSequences() const;
private:
    std::string mInFile;
    std::string mOutFile;
    unsigned int mSeed;
    std::size_t mRepetitions;
    std::vector<dawg::alignment> mAlignments;
    mutt mRng;
    trick mTrickster;
    dawg::matic mKimura;
    std::vector<dawg::ma>  mSegmentModels;

    template <typename VectorType>
    void printVectorContents(VectorType &vec) const;
    void printAlignmentInfo(const dawg::alignment &aln) const;

    template <typename Line, typename File>
    void dawgErrorLog(const std::string &msg, Line l, File f) const;
}; // class Dawg

} // namespace dawg

#endif
