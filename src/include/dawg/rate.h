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

namespace dawg {

class rate_model {
public:
	typedef float rate_type;
	
	template<typename It>
	bool create(const std::string &rname, It first, It last, mutt &m) {
		static std::string name_keys[] = {
			std::string("const"),
			std::string("gamma-invariant")
		};
		switch(key_switch(rname, name_keys)) {
			case 0:
				return create_const(first, last, m);
			case 1:
				return create_gamma(first, last, m);
		};
		return DAWG_ERROR("Invalid rate model; no model named '" << rname << "'");		
		
	}
	
	template<typename It>
	inline bool create_const(It first, It last, mutt &) {
		name = "const";
		std::vector<double> weights(1, 1.0);
		sample.create_inplace(weights);
		values.assign(1, 1.0f);
		return true;
	}

	template<typename It>
	bool create_gamma(It first, It last, mutt &m) {
		if(first == last)
			return DAWG_ERROR("Invalid rate model; gamma-invariant requires at least 1 parameter");
		double alpha = *first++;
		if(alpha < 0.0) 
			return DAWG_ERROR("Invalid rate model; first gamma-invariant parameter '" << alpha << "' is not >= 0.");
		double iota = 0.0;
		if(first != last) {
			iota = *first++;
			if(iota < 0.0 || iota >= 1.0) 
				return DAWG_ERROR("Invalid rate model; second gamma-invariant parameter '" << iota << "' is not [0,1).");
		}
		// construct weights
		std::vector<double> weights(1+DAWG_GAMMA_SAMPLE_SIZE,
			(1.0-iota)/DAWG_GAMMA_SAMPLE_SIZE);
		weights[0] = iota;
		sample.create_inplace(weights);
		
		// construct values
		values.assign(1+DAWG_GAMMA_SAMPLE_SIZE,0.0f);
		double d = 0.0;
		for(std::size_t u = 1; u < values.size(); ++u) {
			values[u] = static_cast<rate_type>(m.rand_gamma(alpha, 1.0/alpha));
			d += values[u];
		}
		// rescale values so that the expected value is exactly 1.0
		d *= (1.0-iota)/DAWG_GAMMA_SAMPLE_SIZE;
		for(std::size_t u = 1; u < values.size(); ++u) {
			values[u] = static_cast<rate_type>(values[u]/d);
		}
		
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
