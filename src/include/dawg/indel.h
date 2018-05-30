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

#include <dawg/utils/specfunc.h>
#include <dawg/utils.h>
#include <dawg/log.h>
#include <dawg/mutt.h>
#include <dawg/utils/aliastable.h>
#include <dawg/error.h>

namespace dawg {

class indel_model {
public:
    typedef std::vector<double> params_type;
    
    indel_model() {
        std::vector<double> a;
        a.push_back(0.0);
        a.push_back(1.0);
        sample.create_inplace(a);
    }
    
    inline double meansize() const { return mean_size; }
    inline double rate() const { return total_rate; }
    inline double upstream_rate() const { return total_upstream_rate; }
    
    // std::numeric_limits<uint32>::max()
    template<typename It1, typename It2, typename It3>
    bool create(It1 first_n, It1 last_n, It2 first_r, It2 last_r,
                It3 first_p, It3 last_p, unsigned int max_size) {
        static std::string name_keys[] = {
            std::string("user"), std::string("geom"),
            std::string("zeta"), std::string("zipf"), std::string("power-law"),
            std::string("yule-simon"), std::string("lavalette")
        };

        if(max_size > std::numeric_limits<boost::uint32_t>::max()) {
	    std::error_code ec = dawg_error::invalid_value;
	    DAWG_ERROR_INFO_ = "maximum indel size is out of range.";
	    throw ec;
        }
        if(first_n == last_n) {
	    std::error_code ec = dawg_error::indel_model_no_type;
	    throw ec;
        }

        total_rate = 0.0;

        for(It2 it=first_r;it!=last_r;++it) {
            if(*it < 0.0) {
		std::error_code ec = dawg_error::invalid_value;
		DAWG_ERROR_INFO_ = "invalid indel model; rate" + boost::lexical_cast<std::string>(*it) + " must be positive.";
		throw ec;
	    }
            total_rate += *it;
        }
        if(total_rate == 0.0) {
            std::vector<double> mix_dist(1, 1.0);
            sample.create(mix_dist);
            sample_upstream.create(mix_dist);
            mean_size = 1.0;
            total_upstream_rate = 0.0;
            return true;
        }

        // enumerate over all models and build table in place
        std::vector<double> mix_dist(max_size+1, 0.0);
        It1 itn = first_n;
        for(It2 it=first_r;it!=last_r;++it) {
            // we have to try to create a model to consume parameters
            bool okay = true;
            double fraction = *it/total_rate;
            switch(key_switch(*itn, name_keys)) {
            case 0: // user model
                okay = create_user(fraction, first_p, last_p, max_size, mix_dist);
                break;
            case 1: // geometric model
                okay = create_geo(fraction, first_p, last_p, max_size, mix_dist);
                break;
            case 2: // zeta power-law model
            case 3:
            case 4: 
                okay = create_zeta(fraction, first_p, last_p, max_size, mix_dist);
                break;
            case 5: // yule-simon model
                okay = create_yule(fraction, first_p, last_p, max_size, mix_dist);
                break;
            case 6: // lavalette model
                okay = create_lavalette(fraction, first_p, last_p, max_size, mix_dist);
                break;
            default:
		std::error_code ec = dawg_error::model_no_name;
		DAWG_ERROR_INFO_ = std::string(itn->c_str()) + " (invalid indel model).";
		throw ec;
            };
            if(!okay) {
		std::error_code ec = dawg_error::creation_fail;
		DAWG_ERROR_INFO_ = "Indel model";
		throw ec;
	    }
            // cycle "name" as needed
            if(++itn == last_n)
                itn = first_n;      
        }
        // Calculate mean
        double m = 0.0, w = 0.0;
        std::vector<double> upstream_dist(mix_dist.size(),0.0);
        for(boost::uint32_t x = static_cast<boost::uint32_t>(mix_dist.size()-1);
                x != 0; --x) {
            m += mix_dist[x]*x;
            w += mix_dist[x];
            upstream_dist[x-1] = upstream_dist[x] + mix_dist[x];
        }
        upstream_dist[0] = 0.0;
        // upstream deletions
        // p(size) = (sum(f(x), size+1, MAX))/(mean(x)-1)
        mean_size = m/w;
        total_upstream_rate = total_rate*(mean_size-1.0);
        sample.create_inplace(mix_dist);
        sample_upstream.create_inplace(upstream_dist);
        
        return true;
    }
        
    boost::uint32_t operator()(mutt &m) const {
        return sample(m.rand_uint64());
    }
    
    boost::uint32_t sample_upstream_overlap(mutt &m) const {
        return sample_upstream(m.rand_uint64());
    }   
private:
    template<typename It>
    inline bool create_geo(double f, It &first, It last, unsigned int max_size,
                           std::vector<double> &mix_dist) {
        if(first == last) {
	    std::error_code ec = dawg_error::param_missing;
	    DAWG_ERROR_INFO_ = "geo requires 1 paramenter (invalid indel model).";
	    throw ec;
	}
        double p = *first++;
        if(p >= 1.0)
            p = 1.0/p;
        if(p <= 0.0) {
	    std::error_code ec = dawg_error::invalid_value;
	    DAWG_ERROR_INFO_ = "geo parameter " + std::to_string(p) + " must be positive (invalid indel model).";
	    throw ec;
	}
        double d = p;
        for(boost::uint32_t n=1; n <= max_size; ++n) {
            mix_dist[n] += f*d;
            d *= (1.0-p);
        }
        return true;
    }

    template<typename It>
    inline bool create_zeta(double f, It &first, It last, unsigned int max_size,
                           std::vector<double> &mix_dist) {
        if(first == last) {
	    std::error_code ec = dawg_error::param_missing;
	    DAWG_ERROR_INFO_ = "zeta requires 1 parameter (invalid indel model).";
	    throw ec;
	}
        double z = *first++;
        if(z <= 1.0) {
	    std::error_code ec = dawg_error::invalid_value;
	    DAWG_ERROR_INFO_ = "zeta parameter " + std::to_string(z) + " must be > 1 (invalid indel model).";
	    throw ec;
	}
        double zz = zeta(z);
        for(boost::uint32_t n=1; n <= max_size; ++n) {
            mix_dist[n] += f*pow(1.0*n,-z)/zz;
        }
        return true;
    }
    
    template<typename It>
    inline bool create_yule(double f, It &first, It last, unsigned int max_size,
                           std::vector<double> &mix_dist) {
        if(first == last) {
	    std::error_code ec = dawg_error::param_missing;
	    DAWG_ERROR_INFO_ = "yule requires 1 parameter (invalid indel model).";
	    throw ec;
	}
        double z = *first++;
        if(z <= 1.0) {
	    std::error_code ec = dawg_error::invalid_value;
	    DAWG_ERROR_INFO_ = "yule parameter " + std::to_string(z) + " must be > 1 (invalid indel model).";
	    throw ec;
	}
        for(boost::uint32_t n=1; n <= max_size; ++n) {
            mix_dist[n] += f*(z-1.0)*beta(n,z);
        }
        return true;
    }

    template<typename It>
    inline bool create_lavalette(double f, It &first, It last, unsigned int max_size,
                           std::vector<double> &mix_dist) {
        if(first == last) {
	    std::error_code ec = dawg_error::param_missing;
	    DAWG_ERROR_INFO_ = "lavalette requires 2 parameters (invalid indel model).";
	    throw ec;
	}
        double z = *first++;
        if(first == last) {
	    std::error_code ec = dawg_error::param_missing;
	    DAWG_ERROR_INFO_ = "lavalette requires 2 parameters (invalid indel model).";
	    throw ec;
	}
        double dm = *first++; 
        if(z <= 1.0) {
	    std::error_code ec = dawg_error::invalid_value;
	    DAWG_ERROR_INFO_ = "lavalette slope " + std::to_string(z) + " must be > 1 (invalid indel model).";
	    throw ec;
	}
        if(dm <= 1.0) {
	    std::error_code ec = dawg_error::invalid_value;
	    DAWG_ERROR_INFO_ = "lavalette max " + std::to_string(dm) + " must be > 1 (invalid indel model).";
	    throw ec;
	}
        boost::uint32_t m = static_cast<boost::uint32_t>(dm);
        // find normalization curve
        double d=0.0;
        for(boost::uint32_t n=m; n != 0; --n) {
            d += pow(1.0*m*n/(m-n+1.0),-z);
        }

        for(boost::uint32_t n=1; n <= m && n <= max_size; ++n) {
            mix_dist[n] += f*pow(1.0*m*n/(m-n+1.0),-z)/d;
        }
        return true;
    }

    
    template<typename It>
    inline bool create_user(double f, It &first, It last, unsigned int max_size,
                           std::vector<double> &mix_dist) {
        double d = 0.0;
        It it = first;
        if(first == last) {
	    std::error_code ec = dawg_error::param_missing;
	    DAWG_ERROR_INFO_ = "no parameteres for user model (invalid indel model).";
	    throw ec;
	}
        // sum up parameters
        for(;it != last && *it >= 0.0;++it)
            d += *it;
        for(boost::uint32_t n=1; first != it && n <= max_size; ++first,++n)
            mix_dist[n] += f*(*first)/d;

        // skip the '-1' terminator if it exists
        if(it != last)
            ++it;
        first = it;
        return true;
    }
    
    alias_table sample, sample_upstream;
    double total_rate, total_upstream_rate;
    double mean_size;
};

} /* namespace dawg */
 
#endif /* DAWG_INDEL_H */
 
