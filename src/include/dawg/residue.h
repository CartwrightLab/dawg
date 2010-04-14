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

namespace dawg {

class residue;

typedef std::vector<residue> sequence;

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
	inline rate_type rate_length() const { return is_deleted() ? rate_type(0.0) : rate_scalar(); }

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
	enum { DNA = 0, RNA, AA, CODON, MODEND, WIDTH_MAX=3 };

	typedef residue_exchange self_type;
	typedef std::pair<const char*,const char*> str_type;

	inline void model(int a, bool lc=false, bool markins=false) {
		if(a >= MODEND)
			return;
		_model = a;
		if(lc) // use lowercase translations
			_model += 3;
		// ??- is a trigraph; add a slash to prevent it
		static const char DNA[] = "ACGT??????????????????????????????????????????????????????????\?-";
		static const char dna[] = "acgt??????????????????????????????????????????????????????????\?-";
		static const char RNA[] = "ACGU??????????????????????????????????????????????????????????\?-";
		static const char rna[] = "acgu??????????????????????????????????????????????????????????\?-";
		// Include O & U, two rare amino acids
		static const char AA[]  = "ACDEFGHIKLMNPQRSTVWYOU????????????????????????????????????????\?-";
		static const char aa[]  = "acdefghiklmnpqrstvwyou????????????????????????????????????????\?-";
		static const char ins[] = "---+++";
		
		static const char* mods[] = { &DNA[0], &RNA[0], &AA[0],
		                              &dna[0], &rna[0], &aa[0]
		                            };
		
		cs_decode = mods[a];
		cs_ins = &ins[(markins ? 3 : 0)];
		width = (a != CODON) ? 1 : 3;
	};

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
			bfirst = std::copy(r.first, r.second, bfirst);
		}
		return bfirst;
	}

	residue_exchange() { model(DNA); }

protected:
	int _model;
	const char* cs_decode;
	const char* cs_ins;
	int width;
};

class residue_factory {
public:
	typedef unsigned int model_type;
	typedef residue_factory self_type;

	static const model_type DNA = 0;
	static const model_type RNA = 1;
	static const model_type AA  = 2;
	static const model_type CODON = 3;
	static const model_type DEL = 0;
	static const model_type DELINS = 1;
	static const model_type INS = 2;

	template<class It, class D>
	void operator()(It b, It e, D &dest) {
		std::transform(b, e, std::back_inserter(dest),
			std::bind1st(std::mem_fun(&self_type::make_residue), this));
	}

	template<class T>
	residue operator()(T a) {
		return make_residue(a);
	}

	struct biop : public std::binary_function<char, double, residue> {
		biop(const residue_factory *o) : obj(o) { }
		const residue_factory *obj;
		residue operator()(char ch, double d) {
			return obj->make_residue(ch,d);
		}
	};
	/*
	template<typename It1, typename It2>
	seq_iterator operator()(It1 b, It1 e, It2 r) {
		store.push_back(Block());
		store.back().reserve(e-b);
		std::transform(b, e, r, std::back_inserter(store.back()), biop(this));
		return store.back().begin();
	}
	*/

	inline void model(model_type a) {_model = a; };

	residue make_residue(char ch) const {
		return residue(encode(ch), static_cast<residue::rate_type>(1.0), 0);
	}
	residue make_residue(char ch, double d) const {
		return residue(encode(ch), static_cast<residue::rate_type>(d), 0);
	}

	residue::data_type encode(char ch) const {
		static residue::data_type dna[] = {0,1,3,2};
		if(_model == DNA || _model == RNA)
			return dna[(ch&6u) >> 1];
		return ~0;
	}

	char decode_gaps(model_type r) const {
		static const char gaps[] = "-=+";
		return gaps[r];
	}

	char decode(residue::data_type r) const {
		static const char dna[] = "ACGT";
		static const char rna[] = "ACGU";
		// Include O & U, two rare amino acids
		static const char aa[] = "ACDEFGHIKLMNPQRSTVWYOU";
		// Encode 64 possible codons using a base-64 notation
		// b/c this function returns a single char.
		// These can later be converted to three nucleotides or 1 amino acid
		static const char cod[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.";

		switch(_model) {
		case DNA:
			return dna[r];
		case RNA:
			return rna[r];
		case AA:
			return aa[r];
		case CODON:
			return cod[r];
		default:
			break;
		}
		return '?';
	}

	residue_factory() : _model(DNA) { }

protected:
	model_type _model;
};

}
#endif

