#pragma once
#ifndef DAWG_RESIDUE_H
#define DAWG_RESIDUE_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>

#include <boost/cstdint.hpp>
#include <boost/bind.hpp>
#include <boost/range/sub_range.hpp>

namespace dawg {

class residue {
public:
	typedef float rate_type;
	typedef boost::uint32_t data_type;

	enum {
		base_mask      =  0x3F, // 00111111
		deleted_base   =  0x3F,
		base_bit_width =  6,
		branch_mask    = ~0x3F,
		branch_inc     =  0x40
	};

	inline data_type base() const { return _data & base_mask; }
	inline void base(data_type b) { _data = (b & base_mask) | (_data & ~base_mask); }

	inline data_type branch() const { return _data & branch_mask; }
	inline void branch(data_type u) { _data = (u & branch_mask) | (_data & ~branch_mask); }

	inline data_type length() const { return is_deleted() ? 0 : 1; }

	inline data_type data()  const { return _data; }
	inline void data(data_type d) { _data = d; }
	inline void data(data_type a, data_type d) {
		_data = (a & base_mask) | (d & branch_mask);
	}

	inline bool is_deleted() const { return (base() == deleted_base); }
	inline bool is_branch(data_type u) const { return (branch() == (u & branch_mask)); }

	inline rate_type rate_scalar() const {return _rate_scalar;}
	inline void rate_scalar(rate_type s) {
		_rate_scalar = s;
	}

	inline data_type rate_cat() const {return _rate_cat;}
	inline void rate_cat(data_type s) {
		_rate_cat = s;
	}

	residue() : _data(0), _rate_scalar(1.0) { }
	residue(data_type xbase, rate_type xscale, data_type xbranch) :
		_data((xbase & base_mask) | (xbranch & branch_mask)),
		_rate_scalar(xscale)
	{

	}

	residue(data_type xbase, data_type xscale, data_type xbranch) :
		_data((xbase & base_mask) | (xbranch & branch_mask)),
		_rate_cat(xscale)
	{

	}


protected:
	data_type  _data;
	union {
	rate_type  _rate_scalar;
	data_type  _rate_cat;
	};
};

template<class CharType, class CharTrait>
std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const dawg::residue &v) {
	if(!o.good()) return o;
	o << v.base();
	return o;
}

class residue_exchange {
public:
	enum { DNA = 0, RNA=2, AA=4, CODON=6, MODEND=30};

	typedef residue_exchange self_type;
	
	typedef boost::sub_range< const char [64] > str_type;

	inline bool model(unsigned int code, bool markins=false, bool keepempty=true) {
		static const char sIns[] = "-+";
		// table for going from base->char
		static const char mods[] =
			"ACGT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // DNA
			"acgt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // dna
			"ACGU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // RNA
			"acgu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // rna
			"ACDEFGHIKLMNPQRSTVWY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // AA
			"acdefghiklmnpqrstvwy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // aa
			"ABCDEFGHIJKLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!!-" // Codons
			"ABCDEFGHIJKLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!!-"
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqruvwxyz0123456789!!!-" 
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-" 
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"ABCDEFGHIJ_:KLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-"
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"ABCDEFGHIJKLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!!-"
			"ABCDEFGHIJKLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!!-"
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"ABCDEFGHIJ_KLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-"
			"ABCDEFGHIJ:KLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"ABCDEFGHIJ:KLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-"
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-"
			"ABCDEFHIJ:KLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!!-"
			"ABDEFGHIJKLNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!!!-"
		;
		// tables for going from char->base
		// 1 genetic code is 80 elements long
		static const char rmods[] = {
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63, 0,63, 1, // DNA & RNA
			63,63,63, 2,63,63,63,63,63,63,63,63,63,63,63,63, 3, 3,63,63,
			63,63,63,63,63,63,63,63,63, 0,63, 1,63,63,63, 2,63,63,63,63,
			63,63,63,63,63,63,63,63, 3, 3,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63, 0,63, 1,
			63,63,63, 2,63,63,63,63,63,63,63,63,63,63,63,63, 3, 3,63,63,
			63,63,63,63,63,63,63,63,63, 0,63, 1,63,63,63, 2,63,63,63,63,
			63,63,63,63,63,63,63,63, 3, 3,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63, 0,63, 1,
			63,63,63, 2,63,63,63,63,63,63,63,63,63,63,63,63, 3, 3,63,63,
			63,63,63,63,63,63,63,63,63, 0,63, 1,63,63,63, 2,63,63,63,63,
			63,63,63,63,63,63,63,63, 3, 3,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63, 0,63, 1,
			63,63,63, 2,63,63,63,63,63,63,63,63,63,63,63,63, 3, 3,63,63,
			63,63,63,63,63,63,63,63,63, 0,63, 1,63,63,63, 2,63,63,63,63,
			63,63,63,63,63,63,63,63, 3, 3,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63, 0,63, 1, // AA
			 2, 3, 4, 5, 6, 7,63, 8, 9,10,11,63,12,13,14,15,16,63,17,18,
			63,19,63,63,63,63,63,63,63, 0,63, 1, 2, 3, 4, 5, 6, 7,63, 8,
			 9,10,11,63,12,13,14,15,16,63,17,18,63,19,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63, 0,63, 1,
			 2, 3, 4, 5, 6, 7,63, 8, 9,10,11,63,12,13,14,15,16,63,17,18,
			63,19,63,63,63,63,63,63,63, 0,63, 1, 2, 3, 4, 5, 6, 7,63, 8,
			 9,10,11,63,12,13,14,15,16,63,17,18,63,19,63,63,63,63,63,63,
			51,52,53,54,55,56,57,58,59,60,63,63,63,63,63,63,63, 0, 1, 2, // Codons
			 3, 4, 5, 6, 7, 8, 9,10,11,63,12,13,14,15,16,17,18,19,20,21,
			22,23,24,63,63,63,63,63,63,25,26,27,28,29,30,31,32,33,34,35,
			36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,63,63,63,63,63,
			51,52,53,54,55,56,57,58,59,60,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,63,12,13,14,15,16,17,18,19,20,21,
			22,23,24,63,63,63,63,63,63,25,26,27,28,29,30,31,32,33,34,35,
			36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,63,63,63,63,63,
			50,51,52,53,54,55,56,57,58,59,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,63,63,44,45,46,47,48,49,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			53,54,55,56,57,58,59,60,61,62,11,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,12,13,63,14,15,16,17,18,19,20,21,22,23,
			24,25,26,63,63,63,63,10,63,27,28,29,30,31,32,33,34,35,36,37,
			38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			51,52,53,54,55,56,57,58,59,60,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,63,12,13,14,15,16,17,18,19,20,21,
			22,23,24,63,63,63,63,63,63,25,26,27,28,29,30,31,32,33,34,35,
			36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,63,63,63,63,63,
			51,52,53,54,55,56,57,58,59,60,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,63,12,13,14,15,16,17,18,19,20,21,
			22,23,24,63,63,63,63,63,63,25,26,27,28,29,30,31,32,33,34,35,
			36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			53,54,55,56,57,58,59,60,61,62,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,11,12,13,14,15,16,17,18,19,20,21,22,23,
			24,25,26,63,63,63,63,10,63,27,28,29,30,31,32,33,34,35,36,37,
			38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,10,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,11,12,63,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,10,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,11,12,63,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,63,
			52,53,54,55,56,57,58,59,60,61,63,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,
			23,24,25,63,63,63,63,63,63,26,27,28,29,30,31,32,33,34,35,36,
			37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,63,63,63,63,63,
			51,52,53,54,55,56,57,58,59,60, 9,63,63,63,63,63,63, 0, 1, 2,
			 3, 4, 5,63, 6, 7, 8,10,11,63,12,13,14,15,16,17,18,19,20,21,
			22,23,24,63,63,63,63,63,63,25,26,27,28,29,30,31,32,33,34,35,
			36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,63,63,63,63,63,
			50,51,52,53,54,55,56,57,58,59,63,63,63,63,63,63,63, 0, 1,63,
			 2, 3, 4, 5, 6, 7, 8, 9,10,63,11,12,13,14,15,16,17,18,19,20,
			21,22,23,63,63,63,63,63,63,24,25,26,27,28,29,30,31,32,33,34,
			35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,63,63,63,63,63
		};
		
		unsigned int a = code%100;
		unsigned int b = code/100;
		
		if(a >= MODEND || mods[a*64] == '!')
			return DAWG_ERROR("invalid genetic code.");
		_model = code;
		_markins = markins;
		_keepempty = keepempty;
				
		cs_decode = &mods[a*64];
		cs_encode = &rmods[a*80];
		cs_ins = &sIns[(_markins ? 1 : 0)];
		
		do_op_append =  (a >= CODON && b < AA) ?
			&residue_exchange::do_op_append_cod :
			&residue_exchange::do_op_append_res ;
		do_op_appendi = (a >= CODON && b < AA) ? 
			&residue_exchange::do_op_appendi_cod :
			&residue_exchange::do_op_appendi_res ;
		
		return true;
	};
	inline bool is_same_model(unsigned int a, bool markins, bool keepempty) const {
		return (a == _model && markins == _markins && keepempty == _keepempty);
	}
	inline bool is_keep_empty() const { return _keepempty; }

	inline residue::data_type encode(char ch) const {
		char ret = (ch >= '0') ? (cs_encode[ch - '0']) : -1;	
		return static_cast<residue::data_type>(ret);
	}
	
	inline char decode(residue::data_type r) const {
		return cs_decode[r & 63];
	}
	
	inline char decode(const residue &r) const {
		return decode(r.base());
	}
	
	inline char decode_ins() const {
		return cs_ins[0];
	}
	
	// codon number -> cod64
	static inline char codon_to_cod64(unsigned int p) {
		const char s[] = "ABCDEFGHIJ_:KLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		return s[p&63];
	}

	// cod64 -> codon number
	static inline unsigned int cod64_to_codon(char c) {
		static const char a[] = {
			// cod64 -> codon number
			54,55,56,57,58,59,60,61,62,63,11,-1,-1,-1,-1,-1,
			-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,12,13,14,15,16,
			17,18,19,20,21,22,23,24,25,26,27,-1,-1,-1,-1,10,
			-1,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
			43,44,45,46,47,48,49,50,51,52,53,-1,-1,-1,-1,-1
		};
		return (c >= '0') ? a[c-'0'] : -1;
	}

	// codon number -> triplet
	static inline unsigned int codon_to_triplet(unsigned int p, int type=0) {
		const char ss[] = "TCAGtcagUCAGucag";
		const char *s = &ss[(type&3)*4];
		unsigned int u = s[p&3];
		u = (u << 8) | s[(p>>2)&3];
		u = (u << 8) | s[(p>>4)&3];
		return u;
	}

	// triplet -> codon number
	static inline unsigned int triplet_to_codon(unsigned int p) {
		// randomly cryptic code:
		// In ncbi descriptions of genetic code,
		//     T|U=0,C=1,A=2,G=3
		// Let f(p) = (p&6)/2.  For ASCII:
		//     f(T|U)= 2, f(C) = 1, f(A)=0, f(G)=3
		// To convert ascii to ncbi: g(p)= (3*p+2)%4
		unsigned int u = p & 0x060606;
		u = (u+(u/2)+0x020202) & 0x030303;
		// reverse bytes and pack into lower 8 bits
		// replaces three shifts, three masks, and two ors
		return ((u*1049601)/65536) % 64;
	}

	void append_residue(std::string &ss, const residue &r) const {
		return (this->*do_op_append)(ss,r);
	}

	void append_ins(std::string &ss) const {
		return (this->*do_op_appendi)(ss);
	}

	static const char* get_protein_code(unsigned int code=0);

	residue_exchange() { model(DNA); }
	

protected:
	void (residue_exchange::*do_op_append)(std::string &ss, const residue &r) const;
	void (residue_exchange::*do_op_appendi)(std::string &ss) const;
	
	void do_op_append_res(std::string &ss, const residue &r) const {
		ss.append(1, decode(r));
	}
	void do_op_append_cod(std::string &ss, const residue &r) const {
		char n = decode(r);
		if(n == '-')
			ss.append(3, '-');
		else {
			unsigned int u = codon_to_triplet(cod64_to_codon(n), _model/100);
			ss.append(1, (char)(u>>16));
			ss.append(1, (char)(u>>8));
			ss.append(1, (char)(u));
		}
	}
	void do_op_appendi_res(std::string &ss) const {
		ss.append(1, decode_ins());
	}
	void do_op_appendi_cod(std::string &ss) const {
		ss.append(3, decode_ins());
	}

	unsigned int _model;
	bool _markins, _keepempty, _translate;
	const char *cs_decode, *cs_ins, *cs_encode;
};

}
#endif

