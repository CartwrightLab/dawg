#pragma once
#ifndef DAWG_DETAILS_SHR3_H
#define DAWG_DETAILS_SHR3_H
/****************************************************************************
 *  Copyright (C) 2009-2012 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#ifndef __STDC_CONSTANT_MACROS
#	define __STDC_CONSTANT_MACROS 1
#endif
#ifndef __STDC_LIMIT_MACROS
#	define __STDC_LIMIT_MACROS 1
#endif
#include <boost/cstdint.hpp>

#if __x86_64__ || __ppc64__ || __amd64__ || __LP64__ || __M_X64 || M_IA64 \
|| _WIN64 || UINTPTR_MAX == 0xffffffffffffffff
#	include <dawg/details/shr3b.h>
#else
#	include <dawg/details/shr3a.h>
#endif
 
#endif
