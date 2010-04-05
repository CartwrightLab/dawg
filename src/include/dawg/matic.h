#pragma once
#ifndef DAWG_MATIC_H
#define DAWG_MATIC_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/ptr_container/ptr_vector.hpp>

#include <dawg/indel.h>
#include <dawg/subst.h>
#include <dawg/rate.h>
#include <dawg/wood.h>
#include <dawg/root.h>
#include <dawg/mutt.h>

#include <dawg/ma.h>

#include <vector>
#include <set>
#include <stack>

namespace dawg {
namespace details {
struct indel_data {
	typedef std::pair<double, boost::uint32_t> element;
	typedef std::stack<element> stack;
	
	stack ins;
	stack del;

	inline void clear() {
		while(!ins.empty())
			ins.pop();
		while(!del.empty())
			del.pop();
	}
};

struct matic_section {
	wood usertree;
	
	subst_model     sub_mod;
	rate_model      rat_mod;
	root_model      rut_mod;
	indel_mix_model ins_mod;
	indel_mix_model del_mod;
	
	void evolve(sequence &child, indel_data &indels, double T, residue::data_type branch_color,
		sequence::const_iterator first, sequence::const_iterator last,
		mutt &m) const;
	void evolve_upstream(sequence &child, indel_data &indels, double T, residue::data_type branch_color,
		sequence::const_iterator first, sequence::const_iterator last,
		mutt &m) const;
	dawg::sequence::const_iterator
	evolve_indels(sequence &child, indel_data &indels, double T, residue::data_type branch_color,
		sequence::const_iterator first, sequence::const_iterator last,
		mutt &m) const;
		
	inline boost::uint32_t next_indel(double d, double &f) const {
		double ins_rate = ins_mod.rate();
		double del_rate = del_mod.rate();
		double indel_rate = ins_rate+del_rate;
		if(d < ins_rate) {
			f = d;
			return 1;
		}
		f = modf((d-ins_rate)/(indel_rate), &d);
		f *= (indel_rate);
		boost::uint32_t x = 2*static_cast<boost::uint32_t>(d);
		if(f < del_rate)
			return 2+x;
		f -= del_rate;
		return 3+x;
	}
};

struct aligned_sequence {
	std::string label;
	std::string seq;
};

}

typedef std::vector<details::aligned_sequence> alignment;

// Core simulation algorithm class
class matic {
public:
	// Configure Simulation
	inline bool configure(const dawg::ma &ma) {
		clear_configuration();
		return add_config_section(ma);
	}
	inline void clear_configuration() {
		configs.clear();
	}
	
	template<class It>
	inline bool configure(It first, It last) {
		clear_configuration();
		for(;first != last; ++first)
			if(!add_config_section(*first))
				return DAWG_ERROR("Configuration section '"
					<< first->name << "' failed to process.");
		return true;
	}

	// Run the simulation
	void walk(alignment& aln);
	
	
	template<typename _It>
	void seed(_It first, _It last) {
		maxx.seed(first, last);
	}
	
	matic() : branch_color(0)
	{}
	
protected:
	typedef dawg::details::matic_section section;
	typedef boost::ptr_vector<section> segment;
	typedef std::vector<segment> segment_vector;
	
	segment_vector configs;
	mutt maxx;
	residue_exchange rex;
	
	residue::data_type branch_color;
	
	bool add_config_section(const dawg::ma &ma);
};

}
#endif //DAWG_MATIC_H
