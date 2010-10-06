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
	enum { DNA = 0, RNA, AA, CODON, MODEND};

	typedef residue_exchange self_type;

	inline void model(int a, int t=0, bool lc=false, bool markins=false, bool keepempty=true, bool translate=false) {
		if(a >= MODEND)
			return;
		_model = a;
		_type = t;
		if(_model == DNA && _type != 1)
			_model = RNA;
		if(lc) // use lowercase translations
			_model += MODEND;
		_markins = markins;
		_keepempty = keepempty;
		_translate = translate;

		static const char DNA[] = "ACGT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-";
		static const char dna[] = "acgt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-";
		static const char RNA[] = "ACGU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-";
		static const char rna[] = "acgu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-";
		static const char AA[]  = "ACDEFGHIKLMNPQRSTVWY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-";
		static const char aa[]  = "acdefghiklmnpqrstvwy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-";
		// codons are complicating depending on translation tables
		// here we will base64 encode codons and translate them somewhere else
		static const char Cod[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-";
		static const char ins[] = "-+";
		//  FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
		// "ABCDEFGHIJ_:KLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
				
		static const char* mods[] = { &DNA[0], &RNA[0], &AA[0], &Cod[0],
		                              &dna[0], &rna[0], &aa[0], &Cod[0]
		                            };
		
		cs_decode = mods[_model];
		cs_ins = &ins[(_markins ? 1 : 0)];
	};
	inline bool is_same_model(int a, int t, bool lc, bool markins, bool keepempty, bool translate) const {
		if(lc)
			a += MODEND;
		return (a == _model && (_model != CODON || t == _type) && markins == _markins
			&& keepempty == _keepempty && translate == _translate);
	}
	inline bool is_keep_empty() const { return _keepempty; }

	inline residue::data_type encode(char ch) const {
		static residue::data_type dna[] = {0,1,3,2};
		static residue::data_type aa[] = {
			20, 0,20, 1, 2, 3, 4, 5, 6, 7,20, 8, 9,10,11,20,
			12,13,14,15,16,20,17,18,20,19,20,20,20,20,20,20
		};
		// char -> triplet num
		static residue::data_type tri[] = {
			54,55,56,57,58,59,60,61,62,63,11,10,10,10,10,10,
			10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,12,13,14,15,16,
			17,18,19,20,21,22,23,24,25,26,27,10,10,10,10,10,
			10,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
			43,44,45,46,47,48,49,50,51,52,53,10,10,10,10,10
		};
		
		switch(_model) {
		case DNA:
		case RNA:
			return dna[(ch&6u) >> 1];
		case AA:
			return aa[((ch&95u)-'@')&31];
		case CODON:
			if('0' <= ch)
				return tri[ ch - '0'];
		default:
		};
			
		return static_cast<residue::data_type>(~0);
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
	
	template<typename It1, typename It2>
	It2 decode_array(It1 afirst, It1 alast, It2 bfirst) const {
		for(;afirst != alast;++afirst) {
			*bfirst++ = decode(*afirst);
		}
		return bfirst;
	}
	
	template<typename It1, typename It2>
	It2 encode_array(It1 afirst, It1 alast, It2 bfirst) const {
	
		for(;afirst != alast;++afirst) {
			*bfirst++ = decode(*afirst);
		}
		return bfirst;
	}
	

	residue_exchange() { model(DNA); }

protected:
	int _model, _type;
	bool _markins, _keepempty, _translate;
	const char* cs_decode;
	const char* cs_ins;
};

}
#endif

