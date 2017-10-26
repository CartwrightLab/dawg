#pragma once
#ifndef DAWG_ROOT_H
#define DAWG_ROOT_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/mutt.h>
#include <dawg/subst.h>
#include <dawg/rate.h>
#include <dawg/sequence.h>

namespace dawg {

class root_model {
public:
	bool create(unsigned int len, const std::string &seq,  const std::vector<double> &rates) {
		root_len = (len == 0) ? seq.size() : len;
		name = "stationary";
		root_seq = seq;
		this->rates = rates;

		if (seq.empty() && rates.empty())
			do_op = &root_model::do_stat;
		else if (seq.empty() && !rates.empty())
			do_op = &root_model::do_user_rates;
		else if (!seq.empty() && rates.empty())
			do_op = &root_model::do_user_seq;
		else
			do_op = &root_model::do_user_seq_rates;

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
		using namespace dawg;
		seq.resize(root_len);
		for(auto i = 0; i != root_len; ++i) {
			seq.at(i).rate_cat(r(m));
			seq.at(i).branch(b);
			if (s.seq_type() == residue_exchange::DNA) {
				residue_exchange rex;
				seq.at(i).base(rex.encode(root_seq.at(i)));
			} else {

			}
		}
	}

	void do_user_rates(sequence &seq, mutt &m, const subst_model &s, const rate_model &r, residue::data_type b) const {
		seq.resize(root_len);
		for(sequence::iterator it=seq.begin(); it != seq.end(); ++it) {
			it->base(s(m));
			it->rate_cat(r(m));
			it->branch(b);
		}
	}

	void do_user_seq_rates(sequence &seq, mutt &m, const subst_model &s, const rate_model &r, residue::data_type b) const {
		seq.resize(root_len);
// TODO add seq calculations
		for(sequence::iterator it=seq.begin(); it != seq.end(); ++it) {
			it->rate_cat(r(m));
			it->branch(b);
		}
	}

	unsigned int root_len;
	std::string name;
	std::string root_seq;
	std::vector<double> rates;
};

} // namespace dawg

#endif // DAWG_ROOT_H
