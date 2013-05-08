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
#include <dawg/utils/aliastable.h>

namespace dawg {

class indel_model {
public:
	typedef std::vector<double> params_type;
	
	indel_model() {
		sample.create(params_type(1,1.0));
	}
	
	inline double meansize() const { return mean;}

	// std::numeric_limits<uint32>::max()
	template<typename It1, typename It2, typename It3>
	bool create(It1 first_n, It1 last_n, It2 first_r, It2 last_r,
	            It3 first_p, It3 last_p, unsigned int max_size) {
		static std::string name_keys[] = {
			std::string("user"), std::string("geom"),
			std::string("zeta"), std::string("zipf"), std::string("power-law")
		};

		if(max_size > std::numeric_limits<uint32>::max())
			return DAWG_ERROR("maximum indel size is out of range");
	    if(first_n == last_n)
	    	return DAWG_ERROR("invalid indel model; no model type specified");

		therate = 0.0;

		for(It2 it=first_r;it!=last_r;++it,++sz) {
			if(*it < 0.0)
				return DAWG_ERROR("invalid indel model; rate '" << *it << "' must be positive");
			therate += *it;
		}
		models.resize(sz);
		
		// enumerate over all models and build table in place
		std::vector<double> mix_dist(max_size, 0);
		It1 itn = first_n;
		for(It2 it=first_r;it!=last_r;++it) {
			// we have to try to create a model to consume parameters
			bool okay = true;
			double fraction = *it/therate;
			switch(key_switch(*itn, name_keys)) {
			case 0: // user model
				okay = create_user(fraction, first, last, max_size, mix_dist);
				break;
			case 1: // geometric model
				okay = create_geo(fraction, first, last, max_size, mix_dist);
				break;
			case 2: // zeta power-law model
			case 3:
			case 4: 
				okay = create_zeta(fraction, first, last, max_size, mix_dist);
				break;
			default:
				return DAWG_ERROR("Invalid indel model; no model named '" << name << "'");
			};
			if(!okay)
				return DAWG_ERROR("indel model creation failed");
			if(++itn == last_n)
				itn = first_n;		
		};
		sample.create_inplace(mix_dist);
		return true;
	}
		
	boost::uint32_t operator()(mutt &m) const {
		return sample(m.rand_uint64());
	}
		
private:
	template<typename It>
	inline bool create_geo(double f, It &first, It last, unsigned int max_size,
	                       std::vector<double> &mix_dist) {
		if(first == last) 
			return DAWG_ERROR("Invalid indel model; geo requires 1 parameter");
		double p = *first++;
		if(p >= 1.0)
			p = 1.0/p;
		if(p <= 0.0)
			return DAWG_ERROR("Invalid indel model; geo parameter '" << p
				<< "' must be positive.");
		double d = p;
		for(unsigned int n=1;n<=max_size;++n) {
			mix_dist[n] += f*d;
			d *= (1.0-p);
		}
		return true;
	}

	template<typename It>
	inline bool create_zeta(double f, It &first, It last, unsigned int max_size,
	                       std::vector<double> &mix_dist) {
		if(first == last) 
			return DAWG_ERROR("Invalid indel model; zeta requires 1 parameter");
		double z = *first++;
		if(z <= 1.0) 
			return DAWG_ERROR("Invalid indel model; zeta parameter '" << qorz
				<< "' must be > 1");
		double d = 1.0;
		double zz = zeta(z);
		for(unsigned int n=1;n<=max_size;++n) {
			mix_dist[n] += f*pow(d,-z)/zz;
			d += 1.0;
		}
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
	
	alias_table sample;
};

} /* namespace dawg */
 
#endif /* DAWG_INDEL_H */
 
