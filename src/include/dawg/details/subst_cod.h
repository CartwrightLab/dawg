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
	
	_model = residue_exchange::CODON + code;
	unsigned int gcode = code%100;
	if(gcode >= 24)
		return DAWG_ERROR("Invalid genetic code.");
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
	remove_stops(gcode, s, freqs);
	
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
	if(first2 != last2) { //+F model
		return create_codgtr("codequ+f", code, s.begin(), s.end(), first2, last2);
	}
	return create_codgtr("codequ", code, s.begin(), s.end(), p.begin(), p.end());
}

template<typename It1, typename It2>
bool subst_model::create_codgy(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> p(64,1.0);
	std::vector<double> s(2016,0.0);
	double kappa, omega;	
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; codgy requires two parameters.");
	kappa = *first1++;
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; codgy requires two parameters.");
	omega = *first1++;
	
	const char *cs_code = residue_exchange::get_protein_code(code);
	
	// create substitution table
	unsigned int u=0;
	for(unsigned int i=0;i<63;++i) {
		for(unsigned int j=i+1;j<64;++j) {
			unsigned int kk = i^j;
			// find number of substitutions
			unsigned int k = ((kk/2)|kk)&21;
			k = (k+k/4+k/16)%4;
			double d = (k>1) ? 0.0 : 1.0;
			// transition or transversion
			if((kk/2+kk/8+kk/32)&1)
				d *= kappa;
			// synonmymous
			if(cs_code[i] != cs_code[j])
				d *= omega;
			s[u++] = d;
		}
	}
	
	if(std::distance(first2,last2) > 4)
		return create_codgtr("codgy+f", code, s.begin(), s.end(), first2, last2);
	double f[4];
	if(first2 == last2)
		f[0] = f[1] = f[2] = f[3] = 0.25;
	else if(!create_freqs("codgy", first2, last2, &f[0], &f[4]))
		return false;
	
	// convert A C G T -> T C A G
	std::swap(f[1], f[3]);
	for(int i=0;i<64;++i)
		p[i] = f[(i)%4]*f[(i/4)%4]*f[(i/16)%4];
	return create_codgtr("codgy", code, s.begin(), s.end(), p.begin(), p.end());
}


} // namesapce dawg

#endif //DAWG_SUBST_COD_H

