#pragma once
#ifndef DAWG_RATE_H
#define DAWG_RATE_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

 #include <string>
 
namespace dawg {

class rate_model {
public:

	template<typename It>
	bool create(const std::string &name, It &first, It last) {
		static std::string name_keys[] = {
			std::string("const"),
			std::string("gamma")
		};
		switch(key_switch(name, name_keys)) {
			case 0:
				return create_const(name, first, last);
			case 1:
				return create_gamma(name, first, last);
		};
		return DAWG_ERROR("Invalid rate model; no model named '" << name << "'");		
		
	}
	template<typename It>
	inline bool create_const(const std::string &name, It &first, It last) {
		name = "const";
		do_op = &rate_model::do_const;
		return true;
	}

	template<typename It>
	bool create_gamma(const std::string &name, It &first, It last) {
		if(first == last)
			return DAWG_ERROR("Invalid rate model; gamma requires at least 1 parameter");
		gamma = *first++;
		if(gamma < 0.0) 
			return DAWG_ERROR("Invalid rate model; first gamma parameter '" << gamma << "' is not >= 0.");
		iota = 0.0;
		if(first != last) {
			iota = *first++;
			if(iota < 0.0 || iota > 1.0) 
				return DAWG_ERROR("Invalid rate model; second gamma parameter '" << iota << "' is not [0,1].");
		}
		if(iota == 0.0) {
			if(gamma == 0.0) {
				name = "const";
				do_op = &rate_model::do_const;
				return true;
			} else if(gamma > 1.0) {
				do_op = &rate_model::do_gamma_low;
			} else {
				do_op = &rate_model::do_gamma_high;
			}
		} else if(gamma == 0.0) {
			do_op = &rate_model::do_iota;
		} else if(gamma > 1.0) {
			do_op = &rate_model::do_iota_gamma_low;
		} else {
			do_op = &rate_model::do_iota_gamma_high;
		}
		name = "gamma";
		return true;
	}
	
	inline const std::string& label() const {
		return name;
	}

	double operator()(mutt &m) const {
		return (this->*do_op)(m);
	}	
	
private:
	// pointer that will hold our method
	double (rate_model::*do_op)(mutt &m) const;

	double do_const(mutt &m) const {
		return 1.0;
	}
	
	double do_iota(mutt &m) const {
		return m.rand_bool(iota) ? 0.0 : 1.0;
	}
	
	double do_gamma_low(mutt &m) const {
		//gamma > 1.0
		return m.rand_gamma_low(1.0/gamma, gamma);
	}
	
	double do_gamma_high(mutt &m) const {
		// gamma <= 1.0
		return m.rand_gamma_high(1.0/gamma, gamma);
	}
	
	double do_iota_gamma_low(mutt &m) const {
		//gamma > 1.0
		return m_rand_bool(iota) ? 0.0 : m.rand_gamma_low(1.0/gamma, gamma);
	}
	
	double do_iota_gamma_high(mutt &m) const {
		// gamma <= 1.0
		return m_rand_bool(iota) ? 0.0 : m.rand_gamma_high(1.0/gamma, gamma);
	}	
	
	std::string name;
	double gamma,iota;
};
 
} /* namespace dawg */
 
#endif /* DAWG_RATE_H */