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
#include <cstring>

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
		base_bit_width =  6,
		branch_mask    = ~0x3F,
		branch_inc     =  0x40
	};

	inline data_type base() const { return _data & base_mask; }
	inline void base(data_type b) { _data = (b & base_mask) | (_data & ~base_mask); }

	inline data_type branch() const { return _data & branch_mask; }
	inline void branch(data_type u) { _data = (u & branch_mask) | (_data & ~branch_mask); }

	inline data_type data()  const { return _data; }
	inline void data(data_type d) { _data = d; }
	inline void data(data_type a, data_type d) {
		_data = (a & base_mask) | (d & branch_mask);
	}

	inline bool is_base(data_type u) const { return (base() == (u & base_mask)); }
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
		// TODO: standardize Root.Code = xyz so that xy picks the type/translation
		// table and z picks upper or lowercase
		// TODO: Allow codons to be translated into aa
		static const char mods[] =
			"ACGT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // DNA
			"acgt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // dna
			"ACGU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // RNA
			"acgu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // rna
			"ACDEFGHIKLMNPQRSTVWY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // AA
			"acdefghiklmnpqrstvwy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // aa
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" // codons
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqr!!uvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ@=KL-MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ@-KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-=KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-=KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEF-HIJ!=KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"AB-DEFGHIJ!!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
		;
		// tables for going from char->base
		// 1 genetic code is 80 elements long
		static const char rmods[] = {
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1, // DNA & RNA
			-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3, 3,-1,-1,
			-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,
			-1,-1,-1,-1,-1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1, // AA
			 2, 3, 4, 5, 6, 7,-1, 8, 9,10,11,-1,12,13,14,15,16,-1,17,18,
			-1,19,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1, 2, 3, 4, 5, 6, 7,-1, 8,
			 9,10,11,-1,12,13,14,15,16,-1,17,18,-1,19,-1,-1,-1,-1,-1,-1,
			54,55,56,57,58,59,60,61,62,63,-1,-1,-1,11,-1,-1,10, 0, 1, 2, // CODON
			 3, 4, 5, 6, 7, 8, 9,12,13,15,16,14,17,18,19,20,21,22,23,24,
			25,26,27,-1,-1,-1,-1,-1,-1,28,29,30,31,32,33,34,35,36,37,38,
			39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,-1,-1,-1,-1,-1
		};
		
		unsigned int a = code%100;
		unsigned int b = code/100;
		
		if(a >= MODEND || mods[a*64] == '!')
			return DAWG_ERROR("invalid genetic code.");
		_model = code;
		_markins = markins;
		_keepempty = keepempty;
				
		cs_decode = &mods[a*64];

		_gap = static_cast<unsigned int>(strchr(cs_decode, '-')-cs_decode);

		if(a < AA)
			cs_encode = &rmods[0*80];
		else if(a < CODON)
			cs_encode = &rmods[1*80];
		else
			cs_encode = &rmods[2*80];

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

	inline residue::data_type gap_base() const { return _gap; }

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
		const char s[] = "ABCDEFGHIJ@=KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		return s[p&63];
	}

	// cod64 -> codon number
	static inline unsigned int cod64_to_codon(char c) {
		static const char a[] = {
			// cod64 -> codon number
			54,55,56,57,58,59,60,61,62,63,-1,-1,-1,11,-1,-1,10, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,12,13,15,16,14,17,18,19,20,21,22,23,24,
			25,26,27,-1,-1,-1,-1,-1,-1,28,29,30,31,32,33,34,35,36,37,38,
			39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,-1,-1,-1,-1,-1
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

	explicit residue_exchange(int m=DNA) { model(m); }
	
	inline static const char* get_protein_code(unsigned int code) {
		static const char s[] = 
			"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
			"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
			"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
		;
		return &s[64*(code%24)];
	}

protected:
	void (residue_exchange::*do_op_append)(std::string &ss, const residue &r) const;
	void (residue_exchange::*do_op_appendi)(std::string &ss) const;
	
	void do_op_append_res(std::string &ss, const residue &r) const {
		ss.append(1, decode(r));
	}
	void do_op_append_cod(std::string &ss, const residue &r) const {
		unsigned int b = r.base()&63;
		if(b == _gap)
			ss.append(3, '-');
		else if(cs_decode[b] == '!')
			ss.append(3, '*');
		else {
			unsigned int u = codon_to_triplet(b, _model/100);
			ss.append(1, (char)(u));
			ss.append(1, (char)(u>>8));
			ss.append(1, (char)(u>>16));
		}
	}
	void do_op_appendi_res(std::string &ss) const {
		ss.append(1, decode_ins());
	}
	void do_op_appendi_cod(std::string &ss) const {
		ss.append(3, decode_ins());
	}

	unsigned int _model, _gap;
	bool _markins, _keepempty, _translate;
	const char *cs_decode, *cs_ins, *cs_encode;
};

}
#endif

