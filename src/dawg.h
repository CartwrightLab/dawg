// dawg.h - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_DAWG_H
#define DAWG_DAWG_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#pragma warning(disable: 4702)

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <map>
#include <functional>

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

// Error Reporting
bool DawgError(const char* csErr, ...);  //always returns false

enum FileFormat { FASTA, NEXUS, PHYLIP, CLUSTAL };
bool SetFormat(FileFormat fmt, int nNum, const char* csBlock);
void DawgIniOutput(std::ostream& os);

const unsigned long FlagOutLowerCase = 1;
const unsigned long FlagOutGapPlus = 2;
const unsigned long FlagOutGapSingleChar = 4;
const unsigned long FlagOutTranslate = 8;

template <class Type> class SumValue
{
private:
	Type m_Sum;
public:
	SumValue () : m_Sum((Type)0) { }
	void operator ( ) ( const Type& elem ) {m_Sum += elem;}
    operator Type() const { return m_Sum; }
};

#ifndef HAVE_GETPID
#	ifdef HAVE__GETPID
#		define getpid _getpid
#	endif
#endif

#endif
