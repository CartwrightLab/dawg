// dawg.h - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_DAWG_H
#define DAWG_DAWG_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#ifndef HAVE_GETPID
#	ifdef HAVE__GETPID
#		define getpid _getpid
#	endif
#endif

#ifdef HAVE_SYS_TYPES_H
#	include <sys/types.h>
#endif

#ifdef HAVE_STDDEF_H
#	include <stddef.h>
#endif

#if !HAVE_MALLOC
extern "C" void *rpl_malloc(size_t n);
#endif

#if !HAVE_REALLOC
extern "C" void *rpl_realloc(void *p, size_t n);
#endif

#ifdef HAVE_TIME_H
#	include <time.h>
#endif
#ifdef HAVE_FLOAT_H
#	include <float.h>
#endif
#ifdef HAVE_STDIO_H
#	include <stdio.h>
#endif
#ifdef HAVE_UNISTD_H
#	include <unistd.h>
#endif
#ifdef HAVE_MATH_H
#	include <math.h>
#endif
#ifdef HAVE_ASSERT_H
#	include <assert.h>
#endif
#ifdef HAVE_STDARG_H
#	include <stdarg.h>
#endif
#ifdef HAVE_PROCESS_H
#	include <process.h>
#endif

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <map>
#include <functional>

// Error Reporting
bool DawgError(const char* csErr, ...);  //always returns false
bool DawgWarn(const char* csErr, ...);  //always returns false

bool SetFormat(unsigned int fmt, int nNum, const char* csBlock);
void DawgIniOutput(std::ostream& os);

// File Formats
const unsigned int FormatFasta = 0;
const unsigned int FormatNexus = 1;
const unsigned int FormatPhylip = 2;
const unsigned int FormatClustal = 3;

// Output Flags
const unsigned int FlagOutLowerCase = 1;
const unsigned int FlagOutGapPlus = 2;
const unsigned int FlagOutGapSingleChar = 4;
const unsigned int FlagOutTranslate = 8;

// Nucleotide Numbers
const int NumAdenine	= 0;
const int NumCytosine	= 1;
const int NumThymine	= 2;
const int NumGuanine	= 3;

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
