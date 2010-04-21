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
		for(int j=i+1;j<20;++j) {
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
	for(int i=0;i<20;++i) {
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
	static const float s[190] = {
		1.02704, 0.738998, 0.0302949, 1.58285, 0.021352, 6.17416, 0.210494, 0.39802,
		0.0467304, 0.0811339, 1.41672, 0.306674, 0.865584, 0.567717, 0.049931, 0.316954,
		0.248972, 0.930676, 0.570025, 0.679371, 0.24941, 0.193335, 0.170135, 0.039437,
		0.127395, 1.05947, 0.0304501, 0.13819, 0.906265, 0.0740339, 0.479855, 2.58443,
		0.088836, 0.373558, 0.890432, 0.323832, 0.397915, 0.384287, 0.0848047, 0.154263,
		2.11517, 0.0613037, 0.499462, 3.17097, 0.257555, 0.893496, 0.390482, 0.103754,
		0.315124, 1.19063, 0.1741, 0.404141, 4.25746, 0.934276, 4.85402, 0.509848,
		0.265256, 5.42942, 0.947198, 0.0961621, 1.12556, 3.95629, 0.554236, 3.01201,
		0.131528, 0.198221, 1.43855, 0.109404, 0.423984, 0.682355, 0.161444, 0.24357,
		0.696198, 0.0999288, 0.556896, 0.415844, 0.171329, 0.195081, 0.908598, 0.0988179,
		0.616783, 5.46947, 0.0999208, 0.330052, 4.29411, 0.113917, 3.8949, 0.869489,
		1.54526, 1.54364, 0.933372, 0.551571, 0.528191, 0.147304, 0.439157, 0.102711,
		0.584665, 2.13715, 0.186979, 5.35142, 0.497671, 0.683162, 0.635346, 0.679489,
		3.0355, 3.37079, 1.40766, 1.07176, 0.704939, 0.545931, 1.34182, 0.740169, 0.31944,
		0.96713, 0.344739, 0.493905, 3.97423, 1.61328, 1.02887, 1.22419, 2.12111, 0.512984,
		0.374866, 0.822765, 0.171903, 0.225833, 0.473307, 1.45816, 1.38698, 0.326622,
		1.51612, 2.03006, 0.795384, 0.857928, 0.554413, 4.37802, 2.00601, 1.00214, 0.152335,
		0.588731, 0.649892, 0.187247, 0.118358, 7.8213, 0.305434, 1.80034, 2.05845, 0.196246,
		0.314887, 0.301281, 0.251849, 0.232739, 1.38823, 0.113133, 0.71707, 0.129767,
		0.156557, 1.52964, 0.336983, 0.262569, 0.212483, 0.137505, 0.665309, 0.515706,
		0.0719167, 0.139405, 0.215737, 1.16392, 0.523742, 0.110864, 0.365369, 0.240735,
		0.543833, 0.325711, 0.196303, 6.45428, 0.103604, 3.87344, 0.42017, 0.133264, 0.398618,
		0.428437, 1.086, 0.216046, 0.22771, 0.381533, 0.786993, 0.291148, 0.31473, 2.48539
	};
	static const float p[20] = {
		0.0866279, 0.0193078, 0.0570451, 0.0580589, 0.0384319, 0.0832518, 0.0244313, 0.048466,
		0.0620286, 0.086209, 0.0195027, 0.0390894, 0.0457631, 0.0367281, 0.043972, 0.0695179,
		0.0610127, 0.0708956, 0.0143859, 0.0352742
	};
	if(first2 != last2) { //+F model
		std::string a(rname);
		a += "+F";
		return create_aagtr(a, &s[0], &s[190], first2, last2);
	}
	return create_aagtr(rname, &s[0], &s[190], &p[0], &p[20]);
}

 
} // namespace dawg
 
#endif // DAWG_SUBST_AA_H
 
