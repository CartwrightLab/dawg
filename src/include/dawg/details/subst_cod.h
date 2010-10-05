#pragma once
#ifndef DAWG_SUBST_COD_H
#define DAWG_SUBST_COD_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/
 
namespace dawg {


// name, followed by params, then freqs
// supports 63 different codons
// the 64th will need to be one unmodeled stop codon
// we use the pos=63 to encode a deletion
template<typename It1, typename It2>
bool subst_model::create_codgtr(const char *mod_name, It1 first1, It1 last1, It2 first2, It2 last2) {
	double d = 0.0;
	int u = 0;
	_model = residue_exchange::CODON;
	// do freqs first
	if(!create_freqs(mod_name, first2, last2, &freqs[0], &freqs[63]))
		return false;
	
	// fill params array
	double params[1953];
	u = 0;
	for(;first1 != last1 && u<1953;++first1,++u) {
		if(*first1 < 0)
			return DAWG_ERROR("Invalid subst model; codgtr parameter #" << u
				<< " '" << *first1 << "' is not >= 0.");
		params[u] = *first1;
	}
	if(u != 1953)
		return DAWG_ERROR("Invalid subst model; codgtr requires 1953 parameters.");
	
	// construct substitution matrix
	// do this locally to enable possible optimizations?
	double s[63][63];
	double rs[63];
	u = 0;
	double aa = 0.0;
	for(int i=0;i<63;++i) {
		s[i][i] = 0.0;
		for(int j=i+1;j<63;++j) {
			s[i][j] = s[j][i] = params[u++];
			aa = std::max(aa,s[i][j]);
		}
	}
	// scale the matrix to substitution time and uniformize
	d = 0.0;
	uni_scale = 0.0;
	for(int i=0;i<63;++i) {
		for(int j=0;j<63;++j) {
			s[i][j] *= freqs[j];
			d += s[i][j]*freqs[i];
		}
	}
	for(int i=0;i<63;++i) {
		rs[i] = 0.0;
		for(int j=0;j<63;++j) {
			s[i][j] /= d;
			rs[i] += s[i][j];
		}
		uni_scale = std::max(uni_scale, rs[i]);
	}
	// create pseudosubstitutions and transition frequencies
	for(int i=0;i<63;++i)
		s[i][i] = uni_scale - rs[i];
	for(int i=0;i<63;++i) {
		for(int j=0;j<63;++j)
			s[i][j] /= uni_scale;
	}
	
	// create cumulative frequencies
	d = 0.0;
	for(int i=0;i<62;++i) {
		d += freqs[i];
		freqs[i] = d;
	}
	// we will include 64 sites in our binary search
	// so fill them with 1.0
	std::fill(&freqs[62],&freqs[64], 1.0);
	for(int i=0;i<63;++i) {
		d = 0.0;
		for(int j=0;j<62;++j) {
			d += s[i][j];
			table[i][j] = d;
		}
		// we will include 32 sites in our binary search
		// so fill them with 1.0
		std::fill(&table[i][62],&table[i][64], 1.0);
	}
	name = mod_name;
	do_op_f = &subst_model::do_codgtr_f;
	do_op_s = &subst_model::do_codgtr_s;
	
	return true;
}

template<typename It1, typename It2>
bool subst_model::create_codequ(const char *, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> s(1953,1.0);
	std::vector<double> p(61,1.0);
	p.push_back(0.0);
	p.push_back(0.0);
	if(first2 != last2) { //+F model
		return create_codgtr("codequ+f", s.begin(), s.end(), first2, last2);
	}
	return create_codgtr("codequ", s.begin(), s.end(), p.begin(), p.end());
}

} // namesapce dawg

#endif //DAWG_SUBST_COD_H

