#pragma once
#ifndef DAWG_SUBST_H
#define DAWG_SUBST_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

 #include <vector>
 
 #include <dawg/mutt.h>
 
 namespace dawg {
 
 class subst_model {
 public:
	typedef uint32_t base_type;
	
	// return random base from stat. dist.
	base_type operator()(mutt &m);
	
	// return random mutant base
	base_type operator()(mutt &m, base_type n);
	
 private:
	// pointer that will hold our method
	base_type (subst_model::*do_op)(mutt &m) const;
 
	// name, followed by params, then freqs
	template<typename It1, typename It2>
	bool create_gtr(const std::string &name, It1 first1, It1 last1, It2 first2, It2 last2) {
		double d = 0.0;
		int u = 0;
		// do freqs first
		for(;first2 != last2 && u<4; ++first2) {
			if(*first2 < 0)
				return DAWG_ERROR("Invalid subst model; gtr frequency #" << u
					<< " '" << *first2 << "' is not >= 0.");
			d += *first2;
			freqs[u] = d;
		}
		// standardize
		for(u=0;u<3++u)
			freqs[u] /= d;
		freqs[3] = 1.0;
		
	}
	
	// must hold at least 64 different characters
	double freqs[64];
	double table[64][64];
 
 } /* namespace dawg */
 
 #endif /* DAWG_SUBST_H */
 