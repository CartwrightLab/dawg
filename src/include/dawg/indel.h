#pragma once
#ifndef DAWG_DAWG_INDEL_H
#define DAWG_DAWG_INDEL_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <boost/cstdint.hpp>

#include <dawg/utils/foreach.h>
#include <dawg/utils/zeta.h>
#include <dawg/utils.h>
#include <dawg/log.h>
#include <dawg/mutt.h>

namespace dawg {

class indel_model {
public:
	typedef std::vector<double> params_type;
	
	indel_model() : do_op(&dawg::indel_model::do_geo), qorz(0.9), mean(10.0), name("geo") { }
	
	inline double meansize() const { return mean;}
	
	template<typename It>
	bool create(const std::string &rname, It &first, It last) {
		static std::string name_keys[] = {
			std::string("user"),
			std::string("geom"),
			std::string("zeta"),
			std::string("zipf"),
			std::string("power-law")
		};
		switch(key_switch(rname, name_keys)) {
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
		
	boost::uint32_t operator()(mutt &m) const {
		return (this->*do_op)(m);
	}
	
	inline const std::string& label() const {
		return name;
	}	
	
private:
	// pointer that will hold our method
	boost::uint32_t (indel_model::*do_op)(mutt &m) const;

	boost::uint32_t do_geo(mutt &m) const {
		return m.rand_geometric_q(qorz);
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
		return (boost::uint32_t)search_binary_cont(udata.begin(), udata.end(), m.rand_uint());
	}
	
	template<typename It>
	inline bool create_geo(It &first, It last) {
		if(first == last) 
			return DAWG_ERROR("Invalid indel model; geo requires 1 parameter");
		qorz = *first++;
		if(qorz >= 1.0)
			qorz = 1.0/qorz;
		if(qorz <= 0.0)
			return DAWG_ERROR("Invalid indel model; geo parameter '" << qorz
				<< "' must be positive.");
		name = "geom";
		mean = 1.0/qorz;
		qorz = -log(1.0-qorz);
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
		if(qorz <= 2.0) //TODO: no quiet; mean is ininite, silently remove upstream stuff
			mean = 1.0;
		else
			mean = zeta(qorz-1.0)/zeta(qorz);
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
		zmax = static_cast<boost::uint32_t>(*first++);
		if(zmax < 1)
			return DAWG_ERROR("Invalid indel model; zeta parameter '" << zmax
				<< "' must be >= 1");
		name = "zipf";
		do_op = &dawg::indel_model::do_zipf;
		// mean is HarmonicNumber[M,z-1]/HarmonicNumber[M,z]
		double a = 0.0,b = 0.0, d=1.0;
		for(boost::uint32_t u=1;u<=zmax;++u,d+=1.0) {
			double t = pow(d,-qorz);
			a += d*t;
			b += t;
		}
		mean = a/b;
		return true;
	}
	
	
	template<typename It>
	inline bool create_user(It &first, It last) {
		udata.assign(1,0);
		double m = 0.0, d=0.0, n=1.0;
		mutt::uint_t mx = std::numeric_limits<mutt::uint_t>::max();
		It it=first;
		// find sum and mean
		for(;it != last && *it >= 0.0;++it) {
			m += *it*n;
			d += *it;
			n += 1.0;
		}
		mean = m/d;
		m = 0.0;
		for(It jt=first;jt != it;++jt) {
			m += *jt;
			udata.push_back(static_cast<mutt::uint_t>(mx*(m/d)));
		}
		if(udata.size() == 1)
			return DAWG_ERROR("Invalid indel model; no parameters for user model.");
		
		udata.back() = mx;
		udata.resize(upper_binary(udata.size()), mx);

		// skip the '-1' terminator if it exists
		if(it != last)
			++it;
		name = "user";
		do_op = &dawg::indel_model::do_user;
		return true;
	}
	
	double qorz, mean;
	std::string name;
	boost::uint32_t zmax;
	std::vector<mutt::uint_t> udata;
};

// TODO: optimize for 1 and 2 components?
class indel_mix_model {
public:
	template<typename It1, typename It2, typename It3>
	bool create(It1 first_n, It1 last_n, It2 first_r, It2 last_r,
	            It3 first_p, It3 last_p) {
	    if(first_n == last_n)
	    	return DAWG_ERROR("invalid indel model; no model type specified");

		double d = 0.0;
		therate = 0.0;
		mean = 0.0;
		std::size_t sz=0,u=0;
		mutt::uint_t mx = std::numeric_limits<mutt::uint_t>::max();

		for(It2 it=first_r;it!=last_r;++it,++sz) {
			if(*it < 0.0)
				return DAWG_ERROR("invalid indel model; rate '" << *it << "' must be positive");
			therate += *it;
		}
		mix.resize(sz, mx);
		models.resize(sz);
		It1 itn = first_n;
		u = 0;
		d = 0.0;
		for(It2 it=first_r;it!=last_r;++it,++u) {
			if(!models[u].create(*itn,first_p, last_p))
				return false;
			if(++itn == last_n)
				itn = first_n;
			mean += (*it/therate)*models[u].meansize();
			d += *it;
			mix[u] = static_cast<mutt::uint_t>((d/therate)*mx);
		}
		u = 0;
		while(u<mix.size()) {
			if(mix[u] > 0.0) {
				++u;
				continue;
			}
			mix.erase(mix.begin()+u);
			models.erase(models.begin()+u);
		}
		if(!mix.empty()) {
			mix.back() = mx;
			mix.resize(upper_binary(sz), mx);
		}
		return true;
	}
	
	boost::uint32_t operator()(mutt &m) const {
		std::size_t x = search_binary_cont(mix.begin(), mix.end(), m.rand_uint());
		return models[x](m);
	}
	double rate() const { return therate; }
	double meansize() const { return mean; }
protected:
	std::vector<mutt::uint_t> mix;
	std::vector<indel_model> models;
	double therate;
	double mean;
};

} /* namespace dawg */
 
#endif /* DAWG_INDEL_H */
 
