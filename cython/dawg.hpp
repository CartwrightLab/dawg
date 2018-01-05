#ifndef DAWG_HPP
#define DAWG_HPP

#include <string>
#include <vector>

#include <dawg/mutt.h>
#include <dawg/trick.h>
#include <dawg/residue.h>

namespace dawg {

class Dawg
{
public:
    explicit Dawg();
    explicit Dawg(const unsigned int s);
    explicit Dawg(const std::string& in,
        const std::string& o,
        const unsigned int r,
        const unsigned int s);
    void run();
    void bark() const;
    unsigned int rand(unsigned int a, unsigned int b);
    void trickStats() const;
private:
    std::string inFile;
    std::string outFile;
    unsigned int reps;
    unsigned int seed;
    std::vector<dawg::alignment> mAlignments;
    mutt rng;
    trick mInput;
}; // class Dawg

} // namespace dawg

#endif
