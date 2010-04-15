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

#include <dawg/utils/foreach.h>

#include <vector>
#include <set>
#include <stack>
#include <map>

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
	typedef std::vector<wood::data_type::size_type> wood_meta_type;

	wood usertree;
	wood_meta_type metatree;
	bool gap_overlap;
	
	subst_model     sub_mod;
	rate_model      rat_mod;
	root_model      rut_mod;
	indel_mix_model ins_mod;
	indel_mix_model del_mod;
	
	void evolve(sequence &child, indel_data &indels, double T, residue::data_type branch_color,
		sequence::const_iterator first, sequence::const_iterator last,
		mutt &m) const;
	void evolve_upstream(sequence &child, indel_data &indels, double T, residue::data_type branch_color,
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

struct sequence_data {
	sequence seq;
	dawg::details::indel_data indels;
};

typedef std::map<std::string,sequence_data> seq_map;

}

class alignment : public std::vector<details::aligned_sequence> {
public:
	
};

template<class CharType, class CharTrait>
inline std::basic_ostream<CharType, CharTrait>&
operator<<(std::basic_ostream<CharType, CharTrait>& o, const alignment &aln) {
	foreach(const alignment::value_type &v, aln) {
		o << v.label << "\t" << v.seq << std::endl;
	}
	return o;
}

// Core simulation algorithm class
class matic {
public:
	// Configure Simulation
	inline bool configure(const dawg::ma &ma) {
		clear_configuration();
		return ( add_config_section(ma) &&
		         finalize_configuration() );
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
		return finalize_configuration();
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
	struct segment : public boost::ptr_vector<section> {
		residue_exchange rex;
	};
	typedef std::vector<segment> segment_vector;
	
	typedef std::map<std::string, wood::data_type::size_type> label_to_index_type;
	typedef std::vector<details::sequence_data> seq_buffers_type;
	
	segment_vector configs;
	mutt maxx;
	
	residue::data_type branch_color;
	
	label_to_index_type label_union;
	
	bool add_config_section(const dawg::ma &ma);
	bool finalize_configuration();
	
	void align(alignment& aln, const seq_buffers_type &seqs, const residue_exchange &rex);
};

}
#endif //DAWG_MATIC_H
