#pragma once
#ifndef DAWG_RATE_H
#define DAWG_RATE_H
/****************************************************************************
 *  Copyright (C) 2009,2013 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

// We will approximate the gamma distribution by taking a large sample.
#define DAWG_GAMMA_SAMPLE_SIZE 4095
 
#include <string>
#include <dawg/utils/aliastable.h>

#include <boost/math/distributions/gamma.hpp>

namespace dawg {

class rate_model {
public:
	typedef float rate_type;
	
	template<typename It>
	bool create(const std::string &rname, It first, It last) {
		static std::string name_keys[] = {
			std::string("const"),
			std::string("gamma-invariant")
		};
		switch(key_switch(rname, name_keys)) {
			case 0:
				return create_const(first, last);
			case 1:
				return create_gamma(first, last);
		};
		return DAWG_ERROR("Invalid rate model; no model named '" << rname << "'");		
		
	}
	
	template<typename It>
	inline bool create_const(It first, It last) {
		name = "const";
		std::vector<double> weights(1, 1.0);
		sample.create_inplace(weights);
		values.assign(1, 1.0f);
		return true;
	}

	template<typename It>
	bool create_gamma(It first, It last) {
		if(first == last)
			return DAWG_ERROR("Invalid rate model; gamma-invariant requires at least 1 parameter");
		double alpha = *first++;
		if(alpha < 0.0)
			return DAWG_ERROR("Invalid rate model; first gamma-invariant parameter '"
				<< alpha << "' is not >= 0.");
		double iota = 0.0;
		if(first != last) {
			iota = *first++;
			if(iota < 0.0 || iota >= 1.0) 
				return DAWG_ERROR("Invalid rate model; second gamma-invariant parameter '"
					<< iota << "' is not [0,1).");
		}
		int sz = DAWG_GAMMA_SAMPLE_SIZE;
		if(first != last) {
			sz = static_cast<int>(*first++);
			// use an upper limit to catch user mistakes.
			if(sz < 1 || sz > 65535)
				return DAWG_ERROR("Invalid rate model; third gamma-invariant parameter '"
					<< sz << "' is not in [1,65535].");
		}
		bool do_median = false;
		if(first != last)
			do_median = (*first++ != 0.0);
		
		// construct weights
		double gw = (1.0-iota)/sz;
		std::vector<double> weights(1+sz, gw);
		weights[0] = iota;
		sample.create_inplace(weights);
		
		// construct values
		values.assign(1+sz,0.0f);
		
		if(do_median) {
			boost::math::gamma_distribution<> gamma_dist(alpha,1.0/alpha);
			for(int k=0; k < sz; ++k)
				values[k+1] = boost::math::quantile(gamma_dist, (2*k+1)/(2.0*sz));
		} else {
			std::vector<double> g(sz);
			boost::math::gamma_distribution<> gamma_dist(alpha,1.0/alpha);
			boost::math::gamma_distribution<> gamma_dist2(alpha+1.0,1.0/alpha);	

			for(int k = 0; k < sz-1; ++k)
				g[k] = boost::math::quantile(gamma_dist, (k+1.0)/sz);
			for(int k = 0; k < sz-1; ++k)
				g[k] = boost::math::cdf(gamma_dist2,g[k]);
			values[1] = g[0]*sz;
			for(int k = 1; k < sz-1; ++k)
				values[k+1] = (g[k]-g[k-1])*sz;
			values[sz] = (1.0-g[sz-2])*sz;
		}
		
		// rescale values so that the expected value is exactly 1.0
		double d = 0.0;
		for(int k=1;k<values.size();++k)
			d += values[k]*gw;
		for(int k=1;k<values.size();++k)
			values[k] /= d;
		
		name = "gamma-invariant";
		return true;
	}
	
	inline const std::string& label() const {
		return name;
	}

	double operator()(mutt &m) const {
		return values[sample(m.rand_uint64())];
	}
	
private:
		
	alias_table sample;
	std::vector<rate_type> values;
	std::string name;
};
 
} /* namespace dawg */
 
#endif /* DAWG_RATE_H */
