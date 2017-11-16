#ifndef DAWG_HPP
#define DAWG_HPP

#include <string>

namespace dawg {

class Dawg
{
public:
    explicit Dawg();
    explicit Dawg(const std::string& in,
        const std::string& o,
        const unsigned int r,
        const unsigned int s);
    void run();
    void bark() const;
private:
    std::string inFile;
    std::string outFile;
    unsigned int reps;
    unsigned int seed;
}; // class Dawg

} // namespace dawg

#endif
