#pragma once
#ifndef DAWG_INDEL_H
#define DAWG_INDEL_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <iostream>
#include <algorithm>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/cstdint.hpp>

#include <dawg/utils.h>
#include <dawg/log.h>
#include <dawg/mutt.h>

#ifndef foreach
#	define foreach BOOST_FOREACH
#endif

namespace dawg {

class indel_model {
public:
	typedef std::vector<double> params_type;
	
	indel_model() : qorz(0.9), name("geo"), do_op(&dawg::indel_model::do_geo) { }
	
	template<typename It>
	bool create(const std::string &name, It &first, It last) {
		static std::string name_keys[] = {
			std::string("user"),
			std::string("geom"),
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
	inline bool create_geo(It &first, It last) {
		if(first == last) 
			return DAWG_ERROR("Invalid indel model; geo requires 1 parameter");
		qorz = 1.0-*first++;
		if(qorz <= 0.0 || qorz >= 1.0)
			return DAWG_ERROR("Invalid indel model; geo parameter '" <<qorz
				<< "' is not between (0,1)");
		name = "geom";
		do_op = &dawg::indel_model::do_geo;
		return true;
	}

	template<typename It>
	inline bool create_zeta(It &first, It last) {
		if(first == last) 
			return DAWG_ERROR("Invalid indel model; zeta requires 1 parameter");
		qorz = *first++;
		if(qorz <= 1.0) 
			return DAWG_ERROR("Invalid indel model; zeta parameter '" << qorz
				<< "' must be > 1");
		name = "zeta";
		do_op = &dawg::indel_model::do_zeta;
		return true;
	}

	template<typename It>
	inline bool create_zipf(It &first, It last) {
		if(first == last || first+1 == last)
			return DAWG_ERROR("Invalid indel model; zeta requires 2 parameters");
		qorz = *first++;
		if(qorz <= 1.0) 
			return DAWG_ERROR("Invalid indel model; zeta parameter '" << qorz
				<< "' must be > 1");
		zmax = static_cast<uint32_t>(*first++);
		if(zmax < 1)
			return DAWG_ERROR("Invalid indel model; zeta parameter '" << zmax
				<< "' must be >= 1");
		name = "zipf";
		do_op = &dawg::indel_model::do_zipf;
		return true;
	}
	
	
	template<typename It>
	inline bool create_user(It &first, It last) {
		udata.assign(1,0.0);
		double d = 0.0;
		while(first != last && *first >= 0.0) {
			d += *first++;
			udata.push_back(d);
		}
		if(*first < 0.0)
			++first;
		if(udata.size() == 1)
			return DAWG_ERROR("Invalid indel model; no parameters for user model.");
		foreach(double &p, udata) {
			p /= d;
		}
		name = "user";
		udata.back() = 1.0;
		do_op = &dawg::indel_model::do_user;
		return true;
	}
	
	boost::uint32_t operator()(mutt &m) const {
		return (this->*do_op)(m);
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
	boost::uint32_t do_user(mutt &m) const {
		return std::distance(udata.begin(), std::lower_bound(
		       udata.begin(), udata.end(), m()));
	}

	std::string name;
	double qorz;
	boost::uint32_t zmax;
	std::vector<double> udata;
};

// TODO: optimize for 1 and 2 components
class indel_mix_model {
public:
	bool create(const std::vector<double> &rates, const std::vector<std::string> &names,
	            const std::vector<double> &params) {
		therate = 0.0;
		mix.reserve(rates.size());
		foreach(double d, rates) {
			if(d <= 0.0)
				return DAWG_ERROR("invalid indel mode; rate '" << d << "' must be positive");
			therate += d;
			mix.push_back(therate);
		}
		std::vector<double> temp(rates);
		models.resize(mix.size());
		std::vector<double>::const_iterator pit = params.begin();
		for(std::size_t u = 0;u<mix.size();++u) {
			mix[u] /= therate;
			if(!models[u].create(names[u % names.size()],pit,params.end()))
				return false;
		}
		mix.back() = 1.0;
		return true;		
	}
	
	boost::uint32_t operator()(mutt &m) const {
		std::size_t x = std::distance(mix.begin(),
		                std::lower_bound(mix.begin(), mix.end(), m()));
		return models[x](m);
	}
	double rate() const { return therate; }
protected:
	std::vector<double> mix;
	std::vector<indel_model> models;
	double therate;
};

} /* namespace dawg */
 
#endif /* DAWG_INDEL_H */
 
