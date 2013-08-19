#pragma once
#ifndef DAWG_SUBST_COD_H
#define DAWG_SUBST_COD_H
/****************************************************************************
 *  Copyright (C) 2009-2012 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/
 
namespace dawg {

// name, followed by params, then freqs
// supports 64 different codons
// stop codons will be silently removed from frequency
template<typename It1, typename It2>
bool subst_model::create_codgtr(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double d = 0.0;
	int u = 0;
	
	type_ = residue_exchange::CODON;
	code_ = code;

	unsigned int gcode = code%100;
	if(gcode >= 24)
		return DAWG_ERROR("Invalid genetic code.");
	// do freqs first
	if(!create_freqs("codgtr", first2, last2, &freqs[0], &freqs[64]))
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
	double rs[64];
	u = 0;
	for(int i=0;i<64;++i) {
		table[i][i] = 0.0;
		for(int j=i+1;j<64;++j) {
			table[i][j] = table[j][i] = params[u++];
		}
	}
	
	// remove stop codons
	remove_stops(gcode, table, freqs);
	
	// scale the matrix to substitution time and uniformize
	d = 0.0;
	uni_scale = 0.0;
	for(int i=0;i<64;++i) {
		for(int j=0;j<64;++j) {
			table[i][j] *= freqs[j];
			d += table[i][j]*freqs[i];
		}
	}
	for(int i=0;i<64;++i) {
		rs[i] = 0.0;
		for(int j=0;j<64;++j) {
			table[i][j] /= d;
			rs[i] += table[i][j];
		}
		uni_scale = std::max(uni_scale, rs[i]);
	}
	// create pseudosubstitutions and transition frequencies
	for(int i=0;i<64;++i)
		table[i][i] = uni_scale - rs[i];
	for(int i=0;i<64;++i) {
		for(int j=0;j<64;++j)
			table[i][j] /= uni_scale;
	}

	if(!create_alias_tables())
		return DAWG_ERROR("unable to create alias tables");
	
	uni_scale *= 3; // adjust for codon mutations
	name = mod_name;
	
	return true;
}

template<typename It1, typename It2>
bool subst_model::create_codgy_equ(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double s[2] = {1.0,1.0};
	return create_codgy("codgy-equ", code, &s[0], &s[2], first2, last2);
}

template<typename It1, typename It2>
bool subst_model::create_codgy(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> p(64,1.0);
	std::vector<double> s(2016,0.0);
	
	std::size_t sz = std::distance(first2, last2);
	std::string mod_namex(mod_name);
	
	if(sz == 4) {
		double df[4];
		mod_namex += "+f4";
		if(!create_freqs(mod_namex.c_str(), first2, last2, &df[0], &df[4]))
			return false;
		// ACGT -> TCAG
		std::swap(df[0], df[2]); std::swap(df[0], df[3]);
		for(int i=0;i<64;++i)
			p[i] = df[(i/16)%4]*df[(i/4)%4]*df[(i)%4];
	} else if(sz == 12) {
		double df[12];
		mod_namex += "+f12";
		if(!create_freqs(mod_namex.c_str(), first2, last2, &df[0], &df[12], 4))
			return false;
		// ACGT -> TCAG
		std::swap(df[0], df[2]);  std::swap(df[0], df[3]);
		// ACGT -> TCAG
		std::swap(df[4], df[6]);  std::swap(df[4], df[7]);
		// ACGT -> TCAG
		std::swap(df[8], df[10]); std::swap(df[8], df[11]);
		for(int i=0;i<64;++i)
			p[i] = df[(i/16)%4]*df[(i/4)%4+4]*df[(i)%4+8];
	} else {
		mod_namex += "+f64";
		if(!create_freqs(mod_namex.c_str(), first2, last2, p.begin(), p.end()))
			return false;
	}
	
	double omega, kappa;
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; " << mod_namex << " requires two parameters.");
	omega = *first1++;
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; " << mod_namex << " requires two parameters.");
	kappa = *first1++;
	
	const char *cs_code = residue_exchange::get_protein_code(code);
	const char *cs_diff = get_codon_diff_upper();

	// create substitution table (T,C,A,G)
	// TC, TA, TG, CA, CG, AG
	double ds[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
	ds[0] = ds[5] = kappa;

	unsigned int u=0;
	for(unsigned int i=0;i<63;++i) {
		for(unsigned int j=i+1;j<64;++j) {
			int d = cs_diff[u];
			s[u] = (d != -1) ? ds[(d%8)] : 0.0;
			if(cs_code[i] != cs_code[j])
				s[u] *= omega;
			++u;
		}
	}
	return create_codgtr(mod_namex.c_str(), code, s.begin(), s.end(), p.begin(), p.end());
}

// name, followed by params, then freqs
// supports 64 different codons
// stop codons will be silently removed from frequency
template<typename It1, typename It2>
bool subst_model::create_codmg_cp(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> p(64,1.0);
	std::vector<double> s(2016,0.0);
	
	std::size_t sz = std::distance(first2, last2);
	std::string mod_namex(mod_name);
	
	// df holds the nucleotide frequences for each position in the format
	// T1, C1, A1, G1, T2, C2, A2, G2, T3, C3, A3, G3
	double df[12];
	
	if(sz == 4) {
		mod_namex += "+f4";
		if(!create_freqs(mod_namex.c_str(), first2, last2, &df[0], &df[4]))
			return false;
		// ACGT -> TCAG
		std::swap(df[0], df[2]); std::swap(df[0], df[3]);
		// fill in the rest
		std::copy(&df[0], &df[4], &df[4]);
		std::copy(&df[0], &df[4], &df[8]);
	} else {
		mod_namex += "+f12";
		if(!create_freqs(mod_namex.c_str(), first2, last2, &df[0], &df[12], 4))
			return false;
		// ACGT -> TCAG
		std::swap(df[0], df[2]);  std::swap(df[0], df[3]);
		// ACGT -> TCAG
		std::swap(df[4], df[6]);  std::swap(df[4], df[7]);
		// ACGT -> TCAG
		std::swap(df[8], df[10]); std::swap(df[8], df[11]);
	}
	
	unsigned int u=0;
	double omega, ds[6], dp[64];
	
	omega = *first1++;
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; codmg-cp requires 71 parameters.");
	for(u=0;first1 != last1 && u < 6;++first1)
		ds[u++] = *first1;
	if(u != 6)
		return DAWG_ERROR("Invalid subst model; codmg-cp requires 71 parameters.");
	// create substitution table (T,C,A,G)
	// TC, TA, TG, CA, CG, AG		
	std::swap(ds[0],ds[4]); std::swap(ds[1],ds[2]);
	std::swap(ds[2],ds[5]); std::swap(ds[3],ds[4]);
	for(u=0;first1 != last1 && u < 64;++first1)
		dp[u++] = *first1;
	if(u != 64)
		return DAWG_ERROR("Invalid subst model; codmg-cp requires 71 parameters.");
		
	// Stationary frequencies
	for(int i=0;i<64;++i)
		p[i] = dp[i]*df[(i/16)%4]*df[(i/4)%4+4]*df[(i)%4+8];
		
	const char *cs_code = residue_exchange::get_protein_code(code);
	const char *cs_diff = get_codon_diff_upper();

	u = 0;
	// i & j represent a codon in 112233 format
	for(unsigned int i=0;i<63;++i) {
		for(unsigned int j=i+1;j<64;++j,++u) {
			int d = cs_diff[u];
			if(d == -1) {
				s[u] = 0.0;
				continue;
			}
			// correct for constant positions
			double f = sqrt(dp[i]*dp[j]);
			switch(d/8) {
			case 0: // First position (high bits) differs
				f *= df[(j/4)%4+4]*df[(j)%4+8];
				break;
			case 1: // Second position (mid bits) differs
				f *= df[(j/16)%4]*df[(j)%4+8];
				break;
			case 2: // Third position (low bits) differs
				f *= df[(j/16)%4]*df[(j/4)%4+4];
				break;
			default:
				return DAWG_ERROR("Unexpected error.");
			}
			if(f > 0.0)
				s[u] = ds[(d%8)]/f;
			else
				s[u] = 1.0;
			if(cs_code[i] != cs_code[j])
				s[u] *= omega;
		}
	}
	return create_codgtr(mod_namex.c_str(), code, s.begin(), s.end(), p.begin(), p.end());
}

template<typename It1, typename It2>
bool subst_model::create_codmg_equ(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> s(71,1.0);
	return create_codmg_cp("codmg-equ", code, s.begin(), s.end(), first2, last2);
}

// requires 7 parameters
template<typename It1, typename It2>
bool subst_model::create_codmg(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> s(71,1.0);
	unsigned int u;
	for(u=0;first1 != last1 && u < 7;++u,++first1)
		s[u] = *first1;
	if(u != 7)
		return DAWG_ERROR("Invalid subst model; codmg requires 7 parameters.");
	return create_codmg_cp("codmg", code, s.begin(), s.end(), first2, last2);	
}

// requires 27 parameters
template<typename It1, typename It2>
bool subst_model::create_codmg_aap(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	std::vector<double> s(71,1.0);
	unsigned int u;
	for(u=0;first1 != last1 && u<7; ++u,++first1)
		s[u] = *first1;
	if(u != 7)
		return DAWG_ERROR("Invalid subst model; codmg-aap requires 27 parameters.");
		
	double dp[20];
	for(u=0;first1 != last1 && u<20; ++u,++first1)
		dp[u] = *first1;
	if(u != 20)
		return DAWG_ERROR("Invalid subst model; codmg-aap requires 27 parameters.");

	const char *cs_code = residue_exchange::get_protein_code(code);
	residue_exchange rex(residue_exchange::AA);
	for(u=0;u<64;++u) {
		char x = cs_code[u];
		s[u+7] = (x == '*') ? 1.0 : dp[rex.encode(x)];
	}
	return create_codmg_cp("codmg-aap", code, s.begin(), s.end(), first2, last2);	
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

/* The upper triangle of a matrix that encodes how the 64 codons differ from
   each other.  Value is -1 if they differ by more than one mutation.  Otherwise
   it is a 5-bit number (aabbb) in which (aa) represents the position of the
   difference, and (bbb) encodes the type of nucleotide difference.
*/
const char* subst_model::get_codon_diff_upper() {
	static const char s[] = {
		16,17,18, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,19,20,-1, 8,-1,-1,-1, 9,-1,
		-1,-1,10,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,21,-1,-1, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 0,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 8,-1,-1,
		-1, 9,-1,-1,-1,10,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,16,17,18,11,-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,19,20,-1,11,-1,-1,-1,
		12,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,21,-1,-1,11,-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,11,-1,-1,-1,12,-1,-1,-1,-1,
		-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,16,17,18,13,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,
		-1,-1,-1,19,20,-1,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1, 2,-1,-1,-1,-1,-1,-1,21,-1,-1,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,13,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,16,17,18,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,19,20,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,21,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,
		16,17,18, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,19,
		20,-1, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,21,-1,-1,
		 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 8,-1,-1,
		-1, 9,-1,-1,-1,10,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,16,17,18,11,-1,-1,-1,12,-1,-1,
		-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,19,20,-1,11,-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1,
		 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,21,-1,-1,11,-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,11,
		-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1,-1,16,17,18,13,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,
		-1,-1,-1,19,20,-1,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,21,-1,-1,13,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		 4,-1,-1,-1,-1,-1,-1,-1,-1,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,16,17,18,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1, 4,-1,-1,-1,19,20,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,21,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1, 4,16,17,18, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 5,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,19,20,-1, 8,-1,-1,-1, 9,-1,
		-1,-1,10,-1,-1,-1, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,21,-1,-1,
		 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1, 5,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,16,17,18,11,-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1, 5,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,19,20,-1,11,-1,-1,-1,12,-1,-1,-1,-1,-1,-1,-1,
		 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,21,-1,-1,11,-1,-1,-1,12,-1,-1,-1,-1,-1,
		-1,-1, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,11,-1,-1,-1,12,-1,-1,-1,-1,
		-1,-1,-1, 5,-1,-1,-1,-1,-1,-1,-1,-1,16,17,18,13,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1, 5,-1,-1,-1,-1,-1,-1,-1,19,20,-1,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1, 5,-1,-1,-1,-1,-1,-1,21,-1,-1,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		 5,-1,-1,-1,-1,-1,-1,-1,-1,13,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,-1,-1,
		-1,-1,16,17,18,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,-1,-1,-1,19,20,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,-1,-1,21,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,
		16,17,18, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1,19,20,-1, 8,-1,-1,-1, 9,-1,
		-1,-1,10,-1,-1,21,-1,-1, 8,-1,-1,-1, 9,-1,-1,-1,10,-1,-1,-1,-1, 8,-1,-1,
		-1, 9,-1,-1,-1,10,16,17,18,11,-1,-1,-1,12,-1,-1,-1,19,20,-1,11,-1,-1,-1,
		12,-1,-1,21,-1,-1,11,-1,-1,-1,12,-1,-1,-1,-1,11,-1,-1,-1,12,16,17,18,13,
		-1,-1,-1,19,20,-1,13,-1,-1,21,-1,-1,13,-1,-1,-1,-1,13,16,17,18,19,20,21
	};
	return &s[0];
}


} // namesapce dawg

#endif //DAWG_SUBST_COD_H

