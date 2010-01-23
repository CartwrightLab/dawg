#pragma once
#ifndef DAWG_ROOT_H
#define DAWG_ROOT_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/mutt.h>
#include <dawg/subst.h>
#include <dawg/rate.h>
#include <dawg/residue.h>

namespace dawg {

class root_model {
public:
	bool create(unsigned int len, const std::string &seq,  const std::vector<double> &rates) {
		root_len = len;
		name = "stat";
		do_op = &root_model::do_stat;
		return true;
	}

	inline void operator()(sequence &seq, mutt &m, const subst_model &s, const rate_model &r) const {
		(this->*do_op)(seq,m,s,r);
	}
	
	inline const std::string& label() const {
		return name;
	}	

private:
	// pointer that will hold our method
	void (root_model::*do_op)(sequence &seq,
		mutt &m, const subst_model &s, const rate_model &r) const;
	
	void do_stat(sequence &seq, mutt &m, const subst_model &s, const rate_model &r) const {
		seq.resize(root_len);
		for(sequence::iterator it=seq.begin(); it != seq.end(); ++it) {
			it->base(s(m));
			it->rate_scalar(r(m));
		}
		
	}
	
	unsigned int root_len;
	std::string name;
};

} // namespace dawg

#endif // DAWG_ROOT_H

