// dawg.h - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_DAWG_H
#define DAWG_DAWG_H

#include <cstdlib>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <cfloat>
#include <cmath>
#include <cassert>
#include <cstdarg>
#include <cstring>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <map>
#include <functional>
#include <ostream>
#include <iterator>

// Error Reporting
bool DawgError(const char* csErr, ...);  //always returns false
bool DawgWarn(const char* csErr, ...);  //always returns false

// Utility functions
#define CERRORR(err_msg) ((std::cerr << "ERROR: " << err_msg << std::endl), false);
#define CERROR(err_msg) (std::cerr << "ERROR: " << err_msg << std::endl);

bool SetFormat(unsigned int fmt, int nNum,
			   const char* csHead, const char* csBefore,
			   const char* csAfter, const char* csTail,
			   bool bSubst);
void DawgIniOutput(std::ostream& os);
void DawgFinOutput(std::ostream& os);

// File Formats
const unsigned int FormatFasta = 0;
const unsigned int FormatNexus = 1;
const unsigned int FormatPhylip = 2;
const unsigned int FormatClustal = 3;

// Output Flags
const unsigned int FlagOutLowerCase     = 1;  // 00001
const unsigned int FlagOutGapPlus       = 2;  // 00010
const unsigned int FlagOutGapSingleChar = 4;  // 00100
const unsigned int FlagOutTranslate     = 8;  // 01000
const unsigned int FlagOutKeepEmpty     = 16; // 10000

// Nucleotide Numbers
const int NumAdenine	= 0;
const int NumCytosine	= 1;
const int NumGuanine	= 2;
const int NumThymine	= 3;

template <class Type> class SumValue
{
private:
	Type m_Sum;
public:
	SumValue () : m_Sum((Type)0) { }
	void operator ( ) ( const Type& elem ) {m_Sum += elem;}
    operator Type() const { return m_Sum; }
};

#endif
