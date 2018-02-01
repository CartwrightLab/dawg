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
        const unsigned int r,
        const unsigned int s);

    explicit Dawg(const std::string& in,
        const std::string& o,
        const unsigned int r,
        const unsigned int s);
    void run();
    void bark() const;
    unsigned int rand(unsigned int a, unsigned int b);
    void trickStats() const;
private:
    std::string mInFile;
    std::string mOutFile;
    std::size_t mRepetitions;
    unsigned int mSeed;
    std::vector<dawg::alignment> mAlignments;
    mutt mRng;
    trick mTrickster;
    // matic_section mMaticSection;


    void printAlignmentInfo(const dawg::alignment &aln) const;
}; // class Dawg

} // namespace dawg

#endif
