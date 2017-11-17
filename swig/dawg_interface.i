%module dawg
%{
// Can put header files here or function declarations like below
#include <string>
extern void run(const std::string &inFile,
    const std::string &outFile,
    const unsigned int reps,
    const unsigned int seed);
%}
