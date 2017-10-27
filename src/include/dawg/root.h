#pragma once
#ifndef DAWG_ROOT_H
#define DAWG_ROOT_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <algorithm>

#include <dawg/mutt.h>
#include <dawg/subst.h>
#include <dawg/rate.h>
#include <dawg/sequence.h>

namespace dawg {

class root_model {
public:
	bool create(unsigned int len, sequence& root_seq) {
		this->root_seq = root_seq;
		root_len = len;

		if (root_seq.empty()) {
			do_op = &root_model::do_stat;
			name = "stationary";
		}
		else if (!root_seq.empty()) {
			do_op = &root_model::do_user_seq;
			name = "user_seq";
		}

		return true;
	}

	inline void operator()(sequence &seq, mutt &m, const subst_model &s, const rate_model &r,
		residue::data_type b) const {
		(this->*do_op)(seq,m,s,r,b);
	}

	inline const std::string& label() const {
		return name;
	}

private:
	// pointer that will hold our method
	void (root_model::*do_op)(sequence &seq, mutt &m,
		const subst_model &s, const rate_model &r, residue::data_type b) const;

	void do_stat(sequence &seq, mutt &m, const subst_model &s, const rate_model &r, residue::data_type b) const {
		seq.resize(root_len);
		for(sequence::iterator it=seq.begin(); it != seq.end(); ++it) {
			it->base(s(m));
			it->rate_cat(r(m));
			it->branch(b);
		}
	}

	void do_user_seq(sequence &seq, mutt &m, const subst_model &s, const rate_model &r, residue::data_type b) const {
		std::copy(root_seq.begin(), root_seq.end(), std::back_inserter(seq));
		// seq.at(i).rate_cat(r(m));
		// seq.at(i).branch(b);
	}

	unsigned int root_len;
	std::string name;
	sequence root_seq;
	std::vector<double> rates;
};

} // namespace dawg

#endif // DAWG_ROOT_H
