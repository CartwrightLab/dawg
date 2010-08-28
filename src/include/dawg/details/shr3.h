#pragma once
#ifndef DAWG_DETAILS_SHR3_H
#define DAWG_DETAILS_SHR3_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#if __x86_64__ || __amd64__ || __M_X64 || M_IA64
#	include <dawg/details/shr3b.h>
namespace dawg { namespace details {
typedef shr3b_mutt_gen shr3_mutt_gen;
}}
#else
#	include <dawg/details/shr3a.h>
namespace dawg { namespace details {
typedef shr3a_mutt_gen shr3_mutt_gen;
}}
#endif
 
#endif
