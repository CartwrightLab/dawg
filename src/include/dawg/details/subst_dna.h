#pragma once
#ifndef DAWG_SUBST_DNA_H
#define DAWG_SUBST_DNA_H
/****************************************************************************
 *  Copyright (C) 2009-2012 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

namespace dawg {

// name, followed by params, then freqs
template<typename It1, typename It2>
bool subst_model::create_gtr(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double d = 0.0;
	int u = 0;

	type_ = residue_exchange::DNA;
	code_ = code;
	
	// do freqs first
	if(!create_freqs(mod_name, first2, last2, &freqs[0], &freqs[4]))
		return false;
	for(int i=4;i<64;++i)
		freqs[i] = 0.0;
	
	// fill params array
	double params[6];
	u = 0;
	for(;first1 != last1 && u<6;++first1,++u) {
		if(*first1 < 0)
			return DAWG_ERROR("Invalid subst model; gtr parameter #" << u
				<< " '" << *first1 << "' is not >= 0.");
		params[u] = *first1;
	}
	if(u != 6)
		return DAWG_ERROR("Invalid subst model; gtr requires six parameters.");
	
	// construct substitution matrix
	double rs[4];
	table[0][0] = table[1][1] = table[2][2] = table[3][3] = 0.0;
	table[0][1] = table[1][0] = params[0]; // A-C
	table[0][2] = table[2][0] = params[1]; // A-G
	table[0][3] = table[3][0] = params[2]; // A-T
	table[1][2] = table[2][1] = params[3]; // C-G
	table[1][3] = table[3][1] = params[4]; // C-T
	table[2][3] = table[3][2] = params[5]; // G-T
	// scale the matrix to substitution time and uniformize
	d = 0.0;
	uni_scale = 0.0;
	for(int i=0;i<4;++i) {
		for(int j=0;j<4;++j) {
			table[i][j] *= freqs[j];
			d += table[i][j]*freqs[i];
		}
	}
	for(int i=0;i<4;++i) {
		rs[i] = 0.0;
		for(int j=0;j<4;++j) {
			table[i][j] /= d;
			rs[i] += table[i][j];
		}
		uni_scale = std::max(uni_scale, rs[i]);
	}
	// create pseudosubstitutions and transition frequencies
	for(int i=0;i<4;++i)
		table[i][i] = uni_scale - rs[i];
	for(int i=0;i<4;++i) {
		for(int j=0;j<4;++j)
			table[i][j] /= uni_scale;
	}
	// fill in the rest of the table matrix
	for(int i=0;i<4;++i)
		for(int j=4;j<64;++j)
			table[i][j] = 0.0;
	for(int i=4;i<64;++i)
		for(int j=0;j<64;++j)
			table[i][j] = 1.0/64.0;
	
	if(!create_alias_tables())
		return DAWG_ERROR("unable to create alias tables");
	name = mod_name;
	return true;
}
	
template<typename It1, typename It2>
bool subst_model::create_jc(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	// equal rates and frequencies
	static const double ones[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
	return create_gtr("jc", code, &ones[0], &ones[6], &ones[0], &ones[4]);
}
template<typename It1, typename It2>
bool subst_model::create_f81(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	// equal rates and frequencies
	static const double ones[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
	return create_gtr("f81", code, &ones[0], &ones[6], first2, last2);
}
template<typename It1, typename It2>
bool subst_model::create_k2p(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	// equal rates and frequencies
	static const double ones[4] = {1.0,1.0,1.0,1.0};
	double p[6], a, b=0.5;  // this default for b means that a=r if b is not specified
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; k2p requires one or two parameters.");
	a = *first1++;
	if(first1 != last1)
		b = *first1;
	p[1] = p[4] = a;
	p[0] = p[2] = p[3] = p[5] = b;
	return create_gtr("k2p", code, &p[0], &p[6], &ones[0], &ones[4]);
}
template<typename It1, typename It2>
bool subst_model::create_tn(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double p[6], f[4], fr, fy, ay, ar, b;
	// read frequencies
	if(!create_freqs(mod_name, first2, last2, &f[0], &f[4]))
		return false;
	fr = f[0]+f[2];
	fy = f[1]+f[3];
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; " << mod_name << " requires two or three parameters.");
	ay = *first1++;
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; " << mod_name << " requires two or three parameters.");
	ar = *first1++;
	if(first1 == last1) {
		// two parameters
		double R = ay, rho=ar;
		ay = (fr*fy*R-f[0]*f[2]-f[1]*f[3])/
		     (2.0*(1.0+R)*(fy*f[0]*f[2]*rho+fr*f[1]*f[3]));
		ar = rho*ay;
		b = 0.5/(fr*fy*(1.0+R));
		ar = ar/fr+b;
		ay = ay/fy+b;
		
	} else {
		// three parameters
		b = *first1;
	}
	p[1] = ar;
	p[4] = ay;
	p[0] = p[2] = p[3] = p[5] = b;	
	return create_gtr(mod_name, code, &p[0], &p[6], &f[0], &f[4]);
}

template<typename It1, typename It2>
bool subst_model::create_tn_f04(const char *mod_name, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double p[6], f[4], fr, fy, ay, ar, b;
	// read frequencies
	if(!create_freqs(mod_name, first2, last2, &f[0], &f[4]))
		return false;
	fr = f[0]+f[2];
	fy = f[1]+f[3];
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; " << mod_name << " requires two or three parameters.");
	ay = *first1++;
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; " << mod_name << " requires two or three parameters.");
	ar = *first1++;
	if(first1 == last1) {
		// two parameters
		double R = ay, rho=ar;
		ay = (fr*fy*R-f[0]*f[2]-f[1]*f[3])/
		     (2.0*(1.0+R)*(fy*f[0]*f[2]*rho+fr*f[1]*f[3]));
		ar = rho*ay;
		b = 0.5/(fr*fy*(1.0+R));
	} else {
		// three parameters
		b = *first1;			
	}
	p[1] = ar/fr+b;
	p[4] = ay/fy+b;
	p[0] = p[2] = p[3] = p[5] = b;	
	return create_gtr(mod_name, code, &p[0], &p[6], &f[0], &f[4]);
}

template<typename It1, typename It2>
bool subst_model::create_f84(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double p[3];
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; f84 requires one or two parameters.");
	double a = *first1++;
	if(first1 == last1) {
		p[0] = a;
		p[1] = 1.0;
		return create_tn_f04("f84", code, &p[0], &p[2], first2, last2);
	}
	double b = *first1;
	p[0] = p[1] = a;
	p[2] = b;
	return create_tn_f04("f84", code, &p[0], &p[3], first2, last2);
}

template<typename It1, typename It2>
bool subst_model::create_hky(const char *, unsigned int code, It1 first1, It1 last1, It2 first2, It2 last2) {
	double p[3];
	if(first1 == last1)
		return DAWG_ERROR("Invalid subst model; hky requires one or two parameters.");
	double a = *first1++;
	if(first1 == last1) {
		p[0] = a;
		p[1] = 1.0;
		return create_tn("hky", code, &p[0], &p[2], first2, last2);
	}
	double b = *first1;
	p[0] = p[1] = a;
	p[2] = b;
	return create_tn("hky", code, &p[0], &p[3], first2, last2);
}

};
 
#endif // DAWG_SUBST_DNA_H

