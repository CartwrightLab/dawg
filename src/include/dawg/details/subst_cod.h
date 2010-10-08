#pragma once
#ifndef DAWG_SUBST_COD_H
#define DAWG_SUBST_COD_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/
 
namespace dawg {

// name, followed by params, then freqs
// supports 64 different codons
// stop codons will be silently removed from the table
template<typename It1, typename It2>
bool subst_model::create_codgtr(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double d = 0.0;
	int u = 0;
	code = (code==0) ? 0 : code-1;
	if(code > 22)
		return DAWG_ERROR("Invalid genetic code.");
	_model = residue_exchange::CODON + code;
	// do freqs first
	if(!create_freqs(mod_name, first2, last2, &freqs[0], &freqs[64]))
		return false;
	
	// fill params array
	double params[2016];
	u = 0;
	for(;first1 != last1 && u<2016;++first1,++u) {
		if(*first1 < 0)
			return DAWG_ERROR("Invalid subst model; codgtr parameter #" << u
				<< " '" << *first1 << "' is not >= 0.");
		params[u] = *first1;
	}
	if(u != 2016)
		return DAWG_ERROR("Invalid subst model; codgtr requires 2016 parameters.");
	
	// construct substitution matrix
	// do this locally to enable possible optimizations?
	double s[64][64];
	double rs[64];
	u = 0;
	for(int i=0;i<64;++i) {
		s[i][i] = 0.0;
		for(int j=i+1;j<64;++j) {
			s[i][j] = s[j][i] = params[u++];
		}
	}
	
	// remove stop codons
	remove_stops(code, freqs, s);
	
	// scale the matrix to substitution time and uniformize
	d = 0.0;
	uni_scale = 0.0;
	for(int i=0;i<64;++i) {
		for(int j=0;j<64;++j) {
			s[i][j] *= freqs[j];
			d += s[i][j]*freqs[i];
		}
	}
	for(int i=0;i<64;++i) {
		rs[i] = 0.0;
		for(int j=0;j<64;++j) {
			s[i][j] /= d;
			rs[i] += s[i][j];
		}
		uni_scale = std::max(uni_scale, rs[i]);
	}
	// create pseudosubstitutions and transition frequencies
	for(int i=0;i<64;++i)
		s[i][i] = uni_scale - rs[i];
	for(int i=0;i<64;++i) {
		for(int j=0;j<64;++j)
			s[i][j] /= uni_scale;
	}
	
	// create cumulative frequencies
	d = 0.0;
	for(int i=0;i<63;++i) {
		d += freqs[i];
		freqs[i] = d;
	}
	freqs[63] = 1.0;
	for(int i=0;i<64;++i) {
		d = 0.0;
		for(int j=0;j<63;++j) {
			d += s[i][j];
			table[i][j] = d;
		}
		// we will include 32 sites in our binary search
		// so fill them with 1.0
		table[i][63] = 1.0;
	}
	
	
	name = mod_name;
	do_op_f = &subst_model::do_codgtr_f;
	do_op_s = &subst_model::do_codgtr_s;
	
	return true;
}

template<typename It1, typename It2>
bool subst_model::create_codequ(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> s(2016,1.0);
	std::vector<double> p(64,1.0);
	p.push_back(0.0);
	p.push_back(0.0);
	if(first2 != last2) { //+F model
		return create_codgtr("codequ+f", code, s.begin(), s.end(), first2, last2);
	}
	return create_codgtr("codequ", code, s.begin(), s.end(), p.begin(), p.end());
}

} // namesapce dawg

#endif //DAWG_SUBST_COD_H

