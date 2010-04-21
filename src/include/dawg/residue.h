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
	enum { DNA = 0, RNA, AA, CODON, MODEND,
	       WIDTH_MAX=3
	     };

	typedef residue_exchange self_type;
	typedef boost::sub_range< const char [64] > str_type;

	inline void model(int a, bool lc=false, bool markins=false, bool keepempty=true) {
		if(a >= MODEND)
			return;
		_model = a;
		if(lc) // use lowercase translations
			_model += MODEND;
		_markins = markins;
		_keepempty = keepempty;
		// ??- is a trigraph; add a slash to prevent it
		static const char DNA[] = "ACGT??????????????????????????????????????????????????????????\?-";
		static const char dna[] = "acgt??????????????????????????????????????????????????????????\?-";
		static const char RNA[] = "ACGU??????????????????????????????????????????????????????????\?-";
		static const char rna[] = "acgu??????????????????????????????????????????????????????????\?-";
		// Include O & U, two rare amino acids
		static const char AA[]  = "ACDEFGHIKLMNPQRSTVWYOU????????????????????????????????????????\?-";
		static const char aa[]  = "acdefghiklmnpqrstvwyou????????????????????????????????????????\?-";
		static const char ins[] = "---+++";
				
		static const char* mods[] = { &DNA[0], &RNA[0], &AA[0], &AA[0],
		                              &dna[0], &rna[0], &aa[0], &aa[0]
		                            };
		
		cs_decode = mods[_model];
		cs_ins = &ins[(_markins ? 3 : 0)];
		width = (a != CODON) ? 1 : 3;
	};
	inline bool is_same_model(int a, bool lc, bool markins, bool keepempty) const {
		if(lc)
			a += MODEND;
		return (a == _model && markins == _markins && keepempty == _keepempty);
	}
	inline bool is_keep_empty() const { return _keepempty; }

	inline residue::data_type encode(char ch) const {
		static residue::data_type dna[] = {0,1,3,2};
		if(_model == DNA || _model == RNA)
			return dna[(ch&6u) >> 1];
		return ~0;
	}
	
	inline str_type decode(residue::data_type r) const {
		r &= 63;
		return std::make_pair(&cs_decode[r],&cs_decode[r+width]);
	}
	
	inline str_type decode(const residue &r) const {
		return decode(r.base());
	}
	
	inline str_type decode_ins() const {
		return std::make_pair(&cs_ins[0], &cs_ins[width]);
	}
	
	template<typename It1, typename It2>
	It2 decode_array(It1 afirst, It1 alast, It2 bfirst) const {
		for(;afirst != alast;++afirst) {
			str_type r = decode(*afirst);
			bfirst = std::copy(r.begin(), r.end(), bfirst);
		}
		return bfirst;
	}

	residue_exchange() { model(DNA); }

protected:
	int _model;
	bool _markins, _keepempty;
	const char* cs_decode;
	const char* cs_ins;
	int width;
};

}
#endif

