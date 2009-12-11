#pragma once
#ifndef DAWG_INDEL_H
#define DAWG_INDEL_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <iostream>
#include <algorithm>
#include <boost/foreach.hpp>

#include <dawg/utils.h>
#include <dawg/mutt.h>

#ifndef foreach
#	define foreach BOOST_FOREACH
#endif

namespace dawg {

class indel_model {
public:
	typedef std::vector<double> params_type;
	
	template<typename It>
	bool create(const std::string &name, It &first, It Last) {
		static std::string name_keys[] = {
			std::string("user"),
			std::string("geo"),
			std::string("zeta"),
			std::string("zipf"), std::string("pl")
		};
		switch(key_switch(name, name_keys)) {
		case 0: // user model
			return create_user(first, last);
		case 1: // geometric model
			return create_geo(first, last);
		case 2: // zeta power-law model
			return create_zeta(first, last);
		case 3:
		case 4: // zipf power-law model
			return create_zipf(first, last);
		};
		return DAWG_ERROR("Invalid indel model; no model named '" << name << "'");
	}
	
	template<typename It>
	inline bool create_geo(It &first, It Last) {
		if(first == last) 
			return DAWG_ERROR("Invalid indel model; geo requires 1 parameter");
		qorz = 1.0-*first++;
		if(qorz <= 0.0 || qorz >= 1.0)
			return DAWG_ERROR("Invalid indel model; geo parameter '" << porz
				<< "' is not between (0,1)");
		name = "geo";
		do_op = &do_geo;
		return true;
	}

	template<typename It>
	inline bool create_zeta(It &first, It Last) {
		if(first == last) 
			return DAWG_ERROR("Invalid indel model; zeta requires 1 parameter");
		qorz = *first++;
		if(qorz <= 1.0) 
			return DAWG_ERROR("Invalid indel model; zeta parameter '" << porz
				<< "' must be > 1");
		name = "zeta";
		do_op = &do_zeta;
		return true;
	}

	template<typename It>
	inline bool create_zipf(It &first, It Last) {
		return false;
	}
	
	
	template<typename It>
	inline bool create_user(It &first, It Last) {
		return false;
	}
	
	boost::unit32_t operator()(mutt &m) {
		return this->*do_op(m);
	}
	
private:
	// pointer that will hold our method
	boost::uint32_t (indel_model::*do_op)(mutt &m) const;

	boost::uint32_t do_geo(mutt &m) const {
		return m.rand_geometric(qorz);
	}
	boost::uint32_t do_zeta(mutt &m) const {
		return m.rand_zeta(qorz);
	}
	boost::uint32_t do_zipf(mutt &m) const {
		boost::uint32_t r;
		do {
			r = m.rand_zeta(qorz);
		} while(r > zmax);
		return r;
	}

	std::string name;
	double qorz;
	boost::uint32_t zmax;
};

// TODO: optimize for 1 and 2 components
class indel_mix_model {
public:
	bool create(const std::vector<double> &rates, const std::vector<double> &names,
	            const std::vector<double> &params) {
		therate = 0.0;
		mix.reserve(rates.size());
		foreach(double d, rates) {
			if(d <= 0.0)
				return DAWG_ERROR("invalid indel mode; rate '" << d << "' must be positive");
			therate += d;
			mix.push_back(d);
		}
		std::vector<double> temp(rates);
		models.resize(mix.size());
		std::vector<double>::const_iterator nit pit = params.begin();
		for(std::size_t u = 0;u<mix.size();++u) {
			mix[u] /= therate;
			models[u].create(names[u % names.size()],pit,params.end());
			// recycle as needed
			if(pit == params.end())
				pit = params.begin();
		}		
	}
	
	boost::uint32_t operator()(mutt &m) const {
		std::size_t x = std::distance(mix.begin(),
		                std::lower_bound(mix.begin(), mix.end(), m()));
		return models[x]();
	}
	double rate() const { return therate; }
protected:
	std::vector<double> mix;
	std::vector<indel_model> models;
	double therate;
};

} /* namespace dawg */
 
#endif /* DAWG_INDEL_H */
 