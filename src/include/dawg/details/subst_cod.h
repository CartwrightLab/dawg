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
	const char *cs_diff = get_codon_diff_upper();

	// create substitution table (T,C,A,G)
	// TC, TA, TG, CA, CG, AG
	double ds[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
	ds[0] = ds[5] = kappa;

	residue_exchange rex;
	rex.model(7);

	unsigned int u=0;
	for(unsigned int i=0;i<63;++i) {
		for(unsigned int j=i+1;j<64;++j) {
			char d = cs_diff[u];
			s[u] = (d <= 5) ? ds[d] : 0.0;
			if(cs_code[i] != cs_code[j])
				s[u] *= omega;
			++u;
		}
	}
	std::size_t freq_len = std::distance(first2,last2);
	if( freq_len == 12) {
		double df[4];
		// 1st position
		if(!create_freqs("codgy", first2, last2, &df[0], &df[4]))
			return false;
		// convert A C G T -> T C A G
		std::swap(df[0], df[2]);
		std::swap(df[0], df[3]);
		for(int i=0;i<64;++i)
			p[i] = df[(i/16)%4];

		// 2nd position
		if(!create_freqs("codgy", first2+4, last2, &df[0], &df[4]))
			return false;
		std::swap(df[0], df[2]);
		std::swap(df[0], df[3]);
		for(int i=0;i<64;++i)
			p[i] *= df[(i/4)%4];

		// 3rd position
		if(!create_freqs("codgy", first2+8, last2, &df[0], &df[4]))
			return false;
		// convert A C G T -> T C A G
		std::swap(df[0], df[2]);
		std::swap(df[0], df[3]);
		for(int i=0;i<64;++i)
			p[i] *= df[(i)%4];
		return create_codgtr("codgy", code, s.begin(), s.end(), p.begin(), p.end());
	} else if(freq_len == 4) {
		double df[4];
		if(!create_freqs("codgy", first2, last2, &df[0], &df[4]))
			return false;
		// convert A C G T -> T C A G
		std::swap(df[0], df[2]);
		std::swap(df[0], df[3]);
		for(int i=0;i<64;++i)
			p[i] = df[(i)%4]*df[(i/4)%4]*df[(i/16)%4];
		return create_codgtr("codgy", code, s.begin(), s.end(), p.begin(), p.end());
	} else if(freq_len == 0) {
		return create_codgtr("codgy", code, s.begin(), s.end(), p.begin(), p.end());
	}
	return create_codgtr("codgy+f", code, s.begin(), s.end(), first2, last2);
}

void subst_model::remove_stops(unsigned int code, double (&s)[64][64], double (&f)[64]) {
	const char *p = residue_exchange::get_protein_code(code);
	if(p[0] == '!')
		return;
	double a = 0.0;
	for(int n=0;n<64;++n) {
		if(p[n] != '*')
			continue;
		a += f[n];
		f[n] = 0.0;
		for(int i=0;i<64;++i)
			s[n][i] = s[i][n] = 1.0;
	}
	a = 1.0-a;
	for(int i=0;i<64;++i)
		f[i] /= a;
}

const char* subst_model::get_codon_diff_upper() {
	static const char s[] = {
		0,1,2,0,9,9,9,1,9,9,9,2,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,
		9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,3,4,9,0,9,9,9,1,9,
		9,9,2,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,
		9,9,2,9,9,9,9,9,9,9,9,9,9,9,9,9,9,5,9,9,0,9,9,9,1,9,9,9,2,9,9,9,0,9,9,9,
		9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,
		9,9,9,9,9,9,9,9,9,0,9,9,9,1,9,9,9,2,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
		9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,9,9,9,9,9,0,1,2,3,9,9,
		9,4,9,9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,
		9,9,9,9,9,2,9,9,9,9,9,9,9,9,9,9,9,3,4,9,3,9,9,9,4,9,9,9,9,9,9,9,0,9,9,9,
		9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,
		9,9,9,5,9,9,3,9,9,9,4,9,9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,
		9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,4,9,9,9,9,
		9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,
		9,9,9,9,9,9,9,9,0,1,2,5,9,9,9,9,9,9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,
		9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,3,4,9,5,9,9,9,9,9,
		9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,
		9,9,2,9,9,9,9,9,9,5,9,9,5,9,9,9,9,9,9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,
		9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,9,5,9,9,9,9,9,9,
		9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
		9,2,9,9,9,9,0,1,2,9,9,9,9,9,9,9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
		9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,9,3,4,9,9,9,9,9,9,9,9,9,9,9,9,9,
		0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,9,9,5,
		9,9,9,9,9,9,9,9,9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,1,9,9,9,9,9,
		9,9,9,9,9,9,9,9,9,9,2,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,0,9,9,9,9,9,9,9,9,
		9,9,9,9,9,9,9,1,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,0,1,2,0,9,9,9,1,9,9,9,2,
		9,9,9,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,3,
		4,9,0,9,9,9,1,9,9,9,2,9,9,9,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,
		9,9,9,9,9,9,9,9,9,5,9,9,0,9,9,9,1,9,9,9,2,9,9,9,3,9,9,9,9,9,9,9,9,9,9,9,
		9,9,9,9,4,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,0,9,9,9,1,9,9,9,2,9,9,9,3,9,9,
		9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,9,9,9,9,9,0,1,2,3,9,9,9,4,9,9,
		9,9,9,9,9,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,9,9,9,9,3,4,9,
		3,9,9,9,4,9,9,9,9,9,9,9,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,
		9,9,9,5,9,9,3,9,9,9,4,9,9,9,9,9,9,9,3,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,
		9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,4,9,9,9,9,9,9,9,3,9,9,9,9,9,9,9,9,9,9,9,9,
		9,9,9,4,9,9,9,9,9,9,9,9,0,1,2,5,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,9,9,9,9,9,
		9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,3,4,9,5,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,9,9,
		9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,5,9,9,5,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,
		9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,3,9,9,
		9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,0,1,2,9,9,9,9,9,9,9,9,9,9,9,9,3,9,9,
		9,9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,9,3,4,9,9,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,
		9,9,9,9,9,9,9,9,9,9,9,9,4,9,9,5,9,9,9,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,9,9,
		9,9,9,9,9,9,9,9,9,9,4,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,9,9,9,9,9,
		9,9,9,9,9,9,9,4,0,1,2,0,9,9,9,1,9,9,9,2,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,9,
		9,9,9,3,4,9,0,9,9,9,1,9,9,9,2,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,9,9,9,5,9,9,
		0,9,9,9,1,9,9,9,2,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,0,9,9,9,1,9,9,
		9,2,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,9,0,1,2,3,9,9,9,4,9,9,9,9,9,9,9,5,9,9,
		9,9,9,9,9,9,9,9,9,3,4,9,3,9,9,9,4,9,9,9,9,9,9,9,5,9,9,9,9,9,9,9,9,9,9,5,
		9,9,3,9,9,9,4,9,9,9,9,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,9,3,9,9,9,4,9,9,9,9,
		9,9,9,5,9,9,9,9,9,9,9,9,0,1,2,5,9,9,9,9,9,9,9,9,9,9,9,5,9,9,9,9,9,9,9,3,
		4,9,5,9,9,9,9,9,9,9,9,9,9,9,5,9,9,9,9,9,9,5,9,9,5,9,9,9,9,9,9,9,9,9,9,9,
		5,9,9,9,9,9,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,5,9,9,9,9,0,1,2,9,9,9,9,9,9,9,
		9,9,9,9,9,5,9,9,9,3,4,9,9,9,9,9,9,9,9,9,9,9,9,9,5,9,9,5,9,9,9,9,9,9,9,9,
		9,9,9,9,9,9,5,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,5,0,1,2,0,9,9,9,1,9,9,9,2,
		9,9,9,3,4,9,0,9,9,9,1,9,9,9,2,9,9,5,9,9,0,9,9,9,1,9,9,9,2,9,9,9,9,0,9,9,
		9,1,9,9,9,2,0,1,2,3,9,9,9,4,9,9,9,3,4,9,3,9,9,9,4,9,9,5,9,9,3,9,9,9,4,9,
		9,9,9,3,9,9,9,4,0,1,2,5,9,9,9,3,4,9,5,9,9,5,9,9,5,9,9,9,9,5,0,1,2,3,4,5
	};
	return &s[0];
}


} // namesapce dawg

#endif //DAWG_SUBST_COD_H

