#pragma once
#ifndef DAWG_SUBST_H
#define DAWG_SUBST_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <vector>
#include <algorithm>

#include <dawg/mutt.h>
#include <dawg/utils.h>
#include <dawg/log.h>
#include <dawg/residue.h>
 
namespace dawg {
 
class subst_model {
public:
	typedef boost::uint32_t base_type;
		
	// return random base from stat. dist.
	inline base_type operator()(mutt &m) const {
		return (this->*do_op_f)(m);
	}
	// return random mutant base
	inline base_type operator()(mutt &m, base_type n) const {
		return (this->*do_op_s)(m,n);
	}

	inline const std::string& label() const { return name; }
	inline double uniform_scale() const { return uni_scale; }
	inline int seq_type() const { return _model; }
	
	template<typename It1, typename It2>
	bool create(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);
	
private:
	// must hold at least 64 different characters
	double freqs[64];
	double table[64][64];
	double uni_scale;
	std::string name;
	unsigned int _model;

	static void remove_stops(unsigned int code, double (&f)[64], double (&s)[64][64]);

	// pointer that will hold our method
	base_type (subst_model::*do_op_f)(mutt &m) const;
	base_type (subst_model::*do_op_s)(mutt &m, base_type n) const;
 
	inline base_type do_gtr_f(mutt &m) const {
		return (base_type)search_binary_cont(&freqs[0], &freqs[4], m());
	}
	inline base_type do_gtr_s(mutt &m, base_type n) const {
		return (base_type)search_binary_cont(&table[n][0], &table[n][4], m());
	}
	inline base_type do_aagtr_f(mutt &m) const {
		return (base_type)search_binary_cont(&freqs[0], &freqs[32], m());
	}
	inline base_type do_aagtr_s(mutt &m, base_type n) const {
		return (base_type)search_binary_cont(&table[n][0], &table[n][32], m());
	}	
	inline base_type do_codgtr_f(mutt &m) const {
		return (base_type)search_binary_cont(&freqs[0], &freqs[64], m());
	}
	inline base_type do_codgtr_s(mutt &m, base_type n) const {
		return (base_type)search_binary_cont(&table[n][0], &table[n][64], m());
	}	

	// DNA Models
	template<typename It1, typename It2>
	bool create_freqs(const char *mod_name, It1 first1, It1 last1, It2 first2, It2 last2) const;
	
	template<typename It1, typename It2>
	bool create_gtr(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);
	
	template<typename It1, typename It2>
	bool create_jc(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_f81(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_k2p(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_tn(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_tn_f04(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_f84(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_hky(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	// Protein Models
	template<typename It1, typename It2>
	bool create_equ(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_aagtr(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_lg(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_wag(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_wagstar(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);
	
	template<typename It1, typename It2>
	bool create_jtt(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_dayhoff(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_molphy(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	// Codon Models	
	template<typename It1, typename It2>
	bool create_codgtr(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_codequ(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2);	
};

template<typename It1, typename It2>
bool subst_model::create(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	static const char name_keys[][16] = {
		"jc",  "gtr", "k2p", "hky", "f84", "f81", "tn", "tn-f04",
		"equ", "aagtr", "lg", "wag", "wagstar", "jtt-dcmut", "dayhoff-dcmut", "molphy",
		"codequ", "codgtr"
	};
	
	static bool (subst_model::*create_ops[])(const char *, unsigned int, It1, It1, It2, It2) = {
		&subst_model::create_jc,  &subst_model::create_gtr, &subst_model::create_k2p,
		&subst_model::create_hky, &subst_model::create_f84, &subst_model::create_f81,
		&subst_model::create_tn,  &subst_model::create_tn_f04,
		&subst_model::create_equ, &subst_model::create_aagtr, &subst_model::create_lg,
		&subst_model::create_wag, &subst_model::create_wagstar, &subst_model::create_jtt,
		&subst_model::create_dayhoff, &subst_model::create_molphy,
		&subst_model::create_codequ, &subst_model::create_codgtr
	};
	std::size_t pos = key_switch(mod_name, name_keys);
	if(pos == (std::size_t)-1)
		return DAWG_ERROR("Invalid subst model; no model named '" << mod_name << "'");
	return (this->*create_ops[pos])(name_keys[pos], code, first1, last1, first2, last2);
}

template<typename It1, typename It2>
bool subst_model::create_freqs(const char *mod_name, It1 first1, It1 last1, It2 first2, It2 last2) const {
	It2 result = first2;
	double d=0.0;
	for(int u=0;first1 != last1 && result != last2;++first1,++result,++u) {
		if(*first1 < 0)
			return DAWG_ERROR("Invalid subst model; " << mod_name << "frequency #" << u
				<< " '" << *first1 << "' is not >= 0.");
		d += *first1;
		*result = *first1;
	}
	if(result != last2)
		return DAWG_ERROR("Invalid subst model; " << mod_name << " requires "
			<< std::distance(first2, last2) << " frequencies.");
	for(;first2 != last2;++first2)
		*first2 /= d;
	return true;
}

} /* namespace dawg */

#include <dawg/details/subst_dna.h>
#include <dawg/details/subst_aa.h>
#include <dawg/details/subst_cod.h>
 
#endif /* DAWG_SUBST_H */
 
