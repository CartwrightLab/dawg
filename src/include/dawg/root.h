#pragma once
#ifndef DAWG_ROOT_H
#define DAWG_ROOT_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <functional>
#include <algorithm>

#include <dawg/mutt.h>
#include <dawg/subst.h>
#include <dawg/rate.h>
#include <dawg/sequence.h>

namespace dawg {

class root_model {
public:
	bool create(unsigned int len, const std::string &seq,  const std::vector<double> &rates) {
		root_len = len;
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
	void calcBaseFromRootSeq(sequence &seq, const subst_model &s) const {
		if (seq.empty())
			return;
		auto triplet = 0;
		std::for_each(root_seq.begin(), root_seq.end(), [&triplet] (const char c)->void {
			triplet |= c;
			triplet <<= 8;
		});

		// DNA
		if (s.seq_type() == 4) {
			seq.at(0).base(residue_exchange::codon_to_triplet(triplet));
		} else {
			seq.at(0).base(residue_exchange::triplet_to_codon(triplet));
		}
	}

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
		seq.resize(root_len);
		calcBaseFromRootSeq(seq, s);

		for(sequence::iterator it=seq.begin(); it != seq.end(); ++it) {
			it->rate_cat(r(m));
			it->branch(b);
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
		calcBaseFromRootSeq(seq, s);

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
