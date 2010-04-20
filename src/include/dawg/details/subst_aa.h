#pragma once
#ifndef DAWG_SUBST_AA_H
#define DAWG_SUBST_AA_H
/****************************************************************************
 *  Copyright (C) 2009-2010 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/
 
namespace dawg {
 
// name, followed by params, then freqs
template<typename It1, typename It2>
bool subst_model::create_aagtr(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
	double d = 0.0;
	int u = 0;
	_model = residue_exchange::AA;
	// do freqs first
	if(!create_freqs(rname, first2, last2, &freqs[0], &freqs[20]))
		return false;
	
	// fill params array
	double params[190];
	u = 0;
	for(;first1 != last1 && u<190;++first1,++u) {
		if(*first1 < 0)
			return DAWG_ERROR("Invalid subst model; aagtr parameter #" << u
				<< " '" << *first1 << "' is not >= 0.");
		params[u] = *first1;
	}
	if(u != 190)
		return DAWG_ERROR("Invalid subst model; aagtr requires 190 parameters.");
	
	// construct substitution matrix
	// do this locally to enable possible optimizations?
	double s[20][20];
	double rs[20];
	u = 0;
	for(int i=0;i<20;++i) {
		s[i][i] = 0.0;
		for(int j=i+1;j<20;++i) {
			s[i][j] = s[j][i] = params[u++];
		}
	}
	// scale the matrix to substitution time and uniformize
	d = 0.0;
	uni_scale = 0.0;
	for(int i=0;i<20;++i) {
		for(int j=0;j<20;++j) {
			s[i][j] *= freqs[j];
			d += s[i][j]*freqs[i];
		}
	}
	for(int i=0;i<20;++i) {
		rs[i] = 0.0;
		for(int j=0;j<20;++j) {
			s[i][j] /= d;
			rs[i] += s[i][j];
		}
		uni_scale = std::max(uni_scale, rs[i]);
	}
	// create pseudosubstitutions and transition frequencies
	for(int i=0;i<20;++i)
		s[i][i] = uni_scale - rs[i];
	for(int i=0;i<20;++i) {
		for(int j=0;j<20;++j)
			s[i][j] /= uni_scale;
	}
	
	// create cumulative frequencies
	d = 0.0;
	for(int i=0;i<19;++i) {
		d += freqs[i];
		freqs[i] = d;
	}
	freqs[19] = 1.0;
	for(int i=0;i<4;++i) {
		d = 0.0;
		for(int j=0;j<19;++j) {
			d += s[i][j];
			table[i][j] = d;
		}
		table[i][19] = 1.0;
	}
	name = rname;
	do_op_f = &subst_model::do_aagtr_f;
	do_op_s = &subst_model::do_aagtr_s;
	
	return true;
}

template<typename It1, typename It2>
bool subst_model::create_wag(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
	static const double s[] = {1.0};
	static const double p[] = {1.0};
}

 
} // namespace dawg
 
#endif // DAWG_SUBST_AA_H
 
