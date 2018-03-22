%module dawg
%include "std_string.i"
%{
void run(const std::string inFile,
    const std::string outFile,
    const unsigned int reps,
    const unsigned int seed);
%}

void run(const std::string inFile,
    const std::string outFile,
    const unsigned int reps,
    const unsigned int seed);
