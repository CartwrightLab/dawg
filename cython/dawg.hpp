#ifndef DAWG_HPP
#define DAWG_HPP

#include <string>

#include <dawg/mutt.h>
#include <dawg/trick.h>

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
    void printSections() const;
    void printSection(const std::string &s);
    void addSection(const std::string &s);
    void setSection(const std::string &key, const std::string &val,
        const std::string &header = std::string("__initial__"));
private:
    std::string inFile;
    std::string outFile;
    unsigned int reps;
    unsigned int seed;
    mutt rng;
    trick input;
}; // class Dawg

} // namespace dawg

#endif
