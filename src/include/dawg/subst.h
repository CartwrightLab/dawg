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
	typedef uint32_t base_type;
	
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
	bool create(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2);
	
private:
	// must hold at least 64 different characters
	double freqs[64];
	double table[64][64];
	double uni_scale;
	std::string name;
	int _model;

	// pointer that will hold our method
	base_type (subst_model::*do_op_f)(mutt &m) const;
	base_type (subst_model::*do_op_s)(mutt &m, base_type n) const;
 
	inline base_type do_gtr_f(mutt &m) const {
		return search_binary_cont(&freqs[0], &freqs[4], m());
	}
	inline base_type do_gtr_s(mutt &m, base_type n) const {
		return search_binary_cont(&table[n][0], &table[n][4], m());
	}
	inline base_type do_aagtr_f(mutt &m) const {
		return search_binary_cont(&freqs[0], &freqs[20], m());
	}
	inline base_type do_aagtr_s(mutt &m, base_type n) const {
		return search_binary_cont(&table[n][0], &table[n][20], m());
	}	

	template<typename It1, typename It2>
	bool create_freqs(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) const;
	
	template<typename It1, typename It2>
	bool create_gtr(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2);
	
	template<typename It1, typename It2>
	bool create_jc(const std::string &, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_f81(const std::string &, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_k2p(const std::string &, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_tn(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_tn_f04(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_f84(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_hky(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2);

	template<typename It1, typename It2>
	bool create_aagtr(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2);

};

template<typename It1, typename It2>
bool subst_model::create(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
	static std::string name_keys[] = {
		std::string("jc"),  std::string("gtr"),
		std::string("k2p"), std::string("hky"),
		std::string("f84"), std::string("f81"), 
		std::string("tn"), std::string("tn-f04")
	};
	switch(key_switch(rname, name_keys)) {
		case 0:
			return create_jc("jc", first1, last1, first2, last2);
		case 1:
			return create_gtr("gtr", first1, last1, first2, last2);
		case 2:
			return create_k2p("k2p", first1, last1, first2, last2);
		case 3:
			return create_hky("hky", first1, last1, first2, last2);
		case 4:
			return create_f84("f84", first1, last1, first2, last2);
		case 5:
			return create_f81("f81", first1, last1, first2, last2);
		case 6:
			return create_tn("tn", first1, last1, first2, last2);
		case 7:
			return create_tn_f04("tn-f04", first1, last1, first2, last2);
			
	};
	return DAWG_ERROR("Invalid subst model; no model named '" << rname << "'");			
}

template<typename It1, typename It2>
bool subst_model::create_freqs(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) const {
	It2 result = first2;
	double d=0.0;
	for(int u=0;first1 != last1 && result != last2;++first1,++result,++u) {
		if(*first1 < 0)
			return DAWG_ERROR("Invalid subst model; " << rname << "frequency #" << u
				<< " '" << *first1 << "' is not >= 0.");
		d += *first1;
		*result = *first1;
	}
	if(result != last2)
		return DAWG_ERROR("Invalid subst model; " << rname << " requires "
			<< std::distance(first2, last2) << " frequencies.");
	for(;first2 != last2;++first2)
		*first2 /= d;
	return true;
}

} /* namespace dawg */

#include <dawg/details/subst_dna.h>
#include <dawg/details/subst_aa.h>
 
#endif /* DAWG_SUBST_H */
 
