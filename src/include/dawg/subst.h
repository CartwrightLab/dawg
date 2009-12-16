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
 
namespace dawg {
 
class subst_model {
public:
	typedef uint32_t base_type;
	
	// return random base from stat. dist.
	inline base_type operator()(mutt &m) {
		return (this->*do_op_f)(m);
	}
	// return random mutant base
	inline base_type operator()(mutt &m, base_type n) {
		return (this->*do_op_s)(m,n);
	}

	inline const std::string& label() const { return name; }
	inline double uniform_scale() const { return uni_scale; }
	
	template<typename It1, typename It2>
	bool create(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
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
	
private:
	// pointer that will hold our method
	base_type (subst_model::*do_op_f)(mutt &m) const;
	base_type (subst_model::*do_op_s)(mutt &m, base_type n) const;
 
	inline base_type do_gtr_f(mutt &m) const {
		return search_binary_cont(&freqs[0], &freqs[4], m());
	}
	inline base_type do_gtr_s(mutt &m, base_type n) const {
		return search_binary_cont(&table[n][0], &table[n][4], m());
	}
	// name, followed by params, then freqs
	template<typename It1, typename It2>
	bool create_gtr(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
		double d = 0.0;
		int u = 0;
		// do freqs first
		if(!create_freqs(rname, first2, last2, &freqs[0], &freqs[4]))
			return false;
		
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
		// do this locally to enable possible optimizations?
		double s[4][4];
		double rs[4];
		s[0][0] = s[1][1] = s[2][2] = s[3][3] = 0.0;
		s[0][1] = s[1][0] = params[0]; // A-C
		s[0][2] = s[2][0] = params[1]; // A-G
		s[0][3] = s[3][0] = params[2]; // A-T
		s[1][2] = s[2][1] = params[3]; // C-G
		s[1][3] = s[3][1] = params[4]; // C-T
		s[2][3] = s[3][2] = params[5]; // G-T
		// scale the matrix to substitution time and uniformize
		d = 0.0;
		uni_scale = 0.0;
		for(int i=0;i<4;++i) {
			for(int j=0;j<4;++j) {
				s[i][j] *= freqs[j];
				d += s[i][j]*freqs[i];
			}
		}
		for(int i=0;i<4;++i) {
			rs[i] = 0.0;
			for(int j=0;j<4;++j) {
				s[i][j] /= d;
				rs[i] += s[i][j];
			}
			uni_scale = std::max(uni_scale, rs[i]);
		}
		// create pseudosubstitutions and transition frequencies
		for(int i=0;i<4;++i)
			s[i][i] = uni_scale - rs[i];
		for(int i=0;i<4;++i) {
			for(int j=0;j<4;++j)
				s[i][j] /= uni_scale;
		}
		
		// create cumulative frequencies
		d = 0.0;
		for(int i=0;i<3;++i) {
			d += freqs[i];
			freqs[i] = d;
		}
		freqs[3] = 1.0;
		for(int i=0;i<4;++i) {
			d = 0.0;
			for(int j=0;j<3;++j) {
				d += s[i][j];
				table[i][j] = d;
			}
			table[i][3] = 1.0;
		}
		name = rname;
		do_op_f = &subst_model::do_gtr_f;
		do_op_s = &subst_model::do_gtr_s;
		
		return true;
	}
	template<typename It1, typename It2>
	bool create_jc(const std::string &, It1 first1, It1 last1, It2 first2, It2 last2) {
		// equal rates and frequencies
		static const double ones[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
		return create_gtr("jc", &ones[0], &ones[6], &ones[0], &ones[4]);
	}
	template<typename It1, typename It2>
	bool create_f81(const std::string &, It1 first1, It1 last1, It2 first2, It2 last2) {
		// equal rates and frequencies
		static const double ones[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
		return create_gtr("f81", &ones[0], &ones[6], first2, last2);
	}
	template<typename It1, typename It2>
	bool create_k2p(const std::string &, It1 first1, It1 last1, It2 first2, It2 last2) {
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
		return create_gtr("k2p", &p[0], &p[6], &ones[0], &ones[4]);
	}
	template<typename It1, typename It2>
	bool create_tn(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
		double p[6], f[4], fr, fy, d, ay, ar, b;
		int u;
		u = 0;
		d = 0.0;
		// read frequencies
		if(!create_freqs(rname, first2, last2, &f[0], &f[4]))
			return false;
		fr = f[0]+f[2];
		fy = f[1]+f[3];
		if(first1 == last1)
			return DAWG_ERROR("Invalid subst model; tn requires two or three parameters.");
		ay = *first1++;
		if(first1 == last1)
			return DAWG_ERROR("Invalid subst model; tn requires two or three parameters.");
		ar = *first1++;
		if(first1 == last1) {
			// two parameters
			double R = ay, rho=ar;
			ay = (fr*fy*R-f[0]*f[2]-f[1]*f[3])/
			     (2.0*(1.0+R)*(fy*f[0]*f[2]*rho+fr*f[1]*f[3]));
			ar = rho*ay;
			b = 0.5/(fr*fy*(1.0+R));
			ar = ar/fr+b;
			ay = ar/fy+b;
			
		} else {
			// three parameters
			b = *first1;			
		}
		p[1] = ar;
		p[4] = ay;
		p[0] = p[2] = p[3] = p[5] = b;	
		return create_gtr(rname, &p[0], &p[6], &f[0], &f[4]);
	}
	
	template<typename It1, typename It2>
	bool create_tn_f04(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
		double p[6], f[4], fr, fy, d, ay, ar, b;
		int u;
		u = 0;
		d = 0.0;
		// read frequencies
		if(!create_freqs(rname, first2, last2, &f[0], &f[4]))
			return false;
		fr = f[0]+f[2];
		fy = f[1]+f[3];
		if(first1 == last1)
			return DAWG_ERROR("Invalid subst model; tn-f04 requires two or three parameters.");
		ay = *first1++;
		if(first1 == last1)
			return DAWG_ERROR("Invalid subst model; tn-f04 requires two or three parameters.");
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
		return create_gtr(rname, &p[0], &p[6], &f[0], &f[4]);
	}
	template<typename It1, typename It2>
	bool create_f84(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
		double p[3];
		if(first1 == last1)
			return DAWG_ERROR("Invalid subst model; f84 requires one or two parameters.");
		double a = *first1++;
		if(first1 == last1) {
			p[0] = a;
			p[1] = 1.0;
			return create_tn_f04("f84", &p[0], &p[2], first2, last2);
		}
		double b = *first1;
		p[0] = p[1] = a;
		p[2] = b;
		return create_tn_f04("f84", &p[0], &p[3], first2, last2);
	}
	template<typename It1, typename It2>
	bool create_hky(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) {
		double p[3];
		if(first1 == last1)
			return DAWG_ERROR("Invalid subst model; hky requires one or two parameters.");
		double a = *first1++;
		if(first1 == last1) {
			p[0] = a;
			p[1] = 1.0;
			return create_tn("hky", &p[0], &p[2], first2, last2);
		}
		double b = *first1;
		p[0] = p[1] = a;
		p[2] = b;
		return create_tn("hky", &p[0], &p[3], first2, last2);
	}
	
	template<typename It1, typename It2>
	bool create_freqs(const std::string &rname, It1 first1, It1 last1, It2 first2, It2 last2) const {
		It2 result = first2;
		double d;
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
		
	// must hold at least 64 different characters
	double freqs[64];
	double table[64][64];
	double uni_scale;
	std::string name;
};
 
} /* namespace dawg */
 
#endif /* DAWG_SUBST_H */
 
