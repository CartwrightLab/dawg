#ifndef DAWG_HPP
#define DAWG_HPP

#include <string>
#include <vector>
#include <map>

#include <dawg/mutt.h>
#include <dawg/trick.h>
#include <dawg/residue.h>
#include <dawg/matic.h>

namespace dawg {

class Dawg
{
public:
    explicit Dawg();
    explicit Dawg(const unsigned int s);

    explicit Dawg(const std::map<std::string, std::vector<std::string>>,
        const std::string& o,
        const unsigned int seed,
        const unsigned int reps);

    explicit Dawg(const std::string& infile,
        const std::string& outfile,
        const unsigned int seed,
        const unsigned int reps);
    void run();
    void printAlignments();
    void bark() const;
    unsigned int rand(unsigned int a, unsigned int b);
    void trickStats() const;
private:
    std::string mInFile;
    std::string mOutFile;
    unsigned int mSeed;
    std::size_t mRepetitions;
    std::vector<dawg::alignment> mAlignments;
    mutt mRng;
    trick mTrickster;
    // matic_section mMaticSection;


    void printAlignmentInfo(const dawg::alignment &aln) const;
}; // class Dawg

} // namespace dawg

#endif
