/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/wood_parse.h>
#include <dawg/matic.h>
#include <dawg/residue.h>
#include <dawg/log.h>

#include <dawg/utils/foreach.h>
#include <dawg/utils/vecio.h>

#include <boost/phoenix/stl/container.hpp>

using namespace dawg;
using namespace std;

template<class It1, class It2>
It2 has_intersection(It1 first1, It1 last1, It2 first2, It2 last2) {
	while(first1 != last1 && first2 != last2) {
		if(*first1 < *first2)
			++first1;
		else if(*first2 < *first1)
			++first2;
		else
			return first2;
	}
	return last2;
}

bool dawg::matic::add_config_section(const dawg::ma &ma) {
	if(ma.tree_model == "user" && ma.tree_tree.empty())
		return true; // Nothing to evolve
	if(ma.root_segment >= configs.size())
		configs.resize(ma.root_segment+1);
	segment &seg = configs[ma.root_segment];
		
	// construction section_info
	std::auto_ptr<section> info(new section);
	
	// construct models
	if(!info->sub_mod.create(ma.subst_model.c_str(), ma.root_code,
		ma.subst_params.begin(), ma.subst_params.end(),
		ma.subst_freqs.begin(), ma.subst_freqs.end()))
		return DAWG_ERROR("substitution model could not be created.");
	
	if(seg.empty()) { // new segment
		if(!seg.rex.model(info->sub_mod.seq_type(), ma.output_markins, ma.output_keepempty))
			return DAWG_ERROR("failed to create sequence type or format object.");
	} else if(!seg.rex.is_same_model(info->sub_mod.seq_type(), ma.output_markins, ma.output_keepempty)) {
		return DAWG_ERROR("the sequence type or format options of a section is different than its segment.");
	}
	info->gap_base = seg.rex.gap_base();
	
	if(!info->rat_mod.create(ma.subst_rate_model, ma.subst_rate_params.begin(),
		ma.subst_rate_params.end(), maxx))
		return DAWG_ERROR("heterogenous rate model could not be created.");
	if(!info->ins_mod.create(ma.indel_model_ins.begin(), ma.indel_model_ins.end(),
		ma.indel_rate_ins.begin(), ma.indel_rate_ins.end(),
		ma.indel_params_ins.begin(), ma.indel_params_ins.end(), ma.indel_max_ins))
		return DAWG_ERROR("insertion model could not be created.");
	if(!info->del_mod.create(ma.indel_model_del.begin(), ma.indel_model_del.end(),
		ma.indel_rate_del.begin(), ma.indel_rate_del.end(),
		ma.indel_params_del.begin(), ma.indel_params_del.end(), ma.indel_max_del))
		return DAWG_ERROR("deletion model could not be created.");
	if(!info->rut_mod.create(ma.root_length, ma.root_seq, ma.root_rates))
		return DAWG_ERROR("root model could not be created.");
	
	// parse tree and find all named descendant nodes
	if(!info->usertree.parse(ma.tree_tree.begin(), ma.tree_tree.end()))
		return DAWG_ERROR("invalid tree; it failed to parse");
		
	// test whether descendents already exist in this segment
	foreach(section &r, seg) {
		std::set<std::string>::const_iterator it = has_intersection(
			r.usertree.desc_labels().begin(), r.usertree.desc_labels().end(),
			info->usertree.desc_labels().begin(), info->usertree.desc_labels().end());
		if(it != info->usertree.desc_labels().end())
			return DAWG_ERROR("invalid tree; descendent '" <<  *it
			               << "' already exists in segment #"
			               << ma.root_segment);
	}
	
	// Allow gap overlap ?
	info->gap_overlap = ma.root_gapoverlap;
	
	// Tree Scale
	info->tree_scale = ma.tree_scale;
	info->usertree.scale(info->tree_scale);
	
	// find location to insert
	segment::iterator it;
	for(it = seg.begin(); it != seg.end()
	    && !info->usertree.has_desc(it->usertree.root_label()); ++it)
		/*noop*/;
	seg.insert(it, info);
	
	return true;
}

bool dawg::matic::finalize_configuration() {
	label_union.clear();
	label_to_index_type::value_type::second_type u = 1;
	
	// remove segments that lack sections and resize vector
	configs.erase(
		remove_if(configs.begin(), configs.end(),
		phoenix::empty(phoenix::arg_names::arg1)),
		configs.end());
	
	// find the union of all labels in the configuration
	foreach(segment &seg, configs) {
		bool has_root = false;
		foreach(section &sec, seg) {
			label_to_index_type::iterator it;
			// try to insert the root label
			it = label_union.insert(label_to_index_type::value_type(sec.usertree.root_label(),0)).first;
			// if this root hasn't been seen before in this segment, mark it
			if(it->second < u) {
				if(has_root)
					return DAWG_ERROR("segment " << u-1 << " has an extra root ( "
						<< sec.usertree.root_label() << " )");
				has_root = true;
			}
			it->second = u;
			// insert and mark all descendant labels
			it = label_union.begin();
			foreach(const string &lab, sec.usertree.desc_labels()) {
				it = label_union.insert(it, label_to_index_type::value_type(lab,0));
				it->second = u;
			}
		}
		++u; // increment segment number
	}
	// set id's for each label
	u = 0;
	aln_size = 0;
	foreach(label_to_index_type::value_type &lab, label_union) {
		lab.second = u++;
		if(lab.first[0] != '{' && lab.first[0] != '~')
			aln_size = u;
	}
	if(aln_size == 0) {
		return DAWG_ERROR("no sequences to align");
	}
	// copy the id's to the meta information for a tree
	foreach(segment &seg, configs) {
		foreach(section &sec, seg) {
			sec.metatree.resize(sec.usertree.data().size());
			for(u=0;u<sec.usertree.data().size();++u) {
				sec.metatree[u] = label_union[sec.usertree.data()[u].label];
			}
		}
	}
	
	return true;
}

void dawg::matic::pre_walk(alignment& aln) {
	if(configs.empty())
		return;
	aln.resize(aln_size);
	seqs.resize(label_union.size());
	aln.max_label_width = 0;
	aln.max_label_width_14 = 14;
	label_to_index_type::const_iterator it = label_union.begin();
	for(alignment::size_type uu = 0; uu < aln_size;++uu) {
		const string &lab = it->first;
		aln.max_label_width = std::max(aln.max_label_width, lab.length());
		aln.max_label_width_14 = std::max(aln.max_label_width_14, lab.length());
		aln[uu].label = lab;
		++it;
	}
	// TODO: We need to either initialize "empty" segments or skip over them
	// Or prune them
	aln.seq_type = configs[0][0].sub_mod.seq_type();
}

void dawg::matic::walk(alignment& aln) {
	if(configs.empty())
		return;
	// clear alignment
	foreach(alignment::value_type &a, aln) {
		a.seq.clear();
	}

	if(configs[0][0].gap_overlap) {
		// take first seg and do upstream indels
		// root.gapoverlap = false prevents upstream indel creation		
		branch_color = 0;
		// clear sequence buffer
		foreach(details::sequence_data &v, seqs) {
			v.seq.clear();
			v.indels.clear();
		}
		foreach(const section &sec, configs[0]) {
			for(wood::data_type::size_type u=1;u<sec.usertree.data().size();++u) {
				const wood::node &n = sec.usertree.data()[u];
				const sequence &ranc = seqs[sec.metatree[u-n.anc]].seq;
				details::sequence_data &sd = seqs[sec.metatree[u]];
				sec.evolve_upstream(sd.seq, sd.indels, n.length,
					branch_color, maxx);
				sec.evolve(sd.seq, sd.indels, n.length,
					branch_color, ranc.begin(), ranc.end(), maxx);
				branch_color += dawg::residue::branch_inc;
			}
		}
		align(aln, seqs, configs[0].rex);
	}

	foreach(const segment &seg, configs) {
		if(seg.empty())
			continue;
		branch_color = 0;
		// clear the sequence buffers
		foreach(details::sequence_data &v, seqs) {
			v.seq.clear();
		}
		{ // create root sequence of this 
			const section &sec = seg[0];
			seqs[sec.metatree[0]].indels.clear();	
			sec.rut_mod(seqs[sec.metatree[0]].seq, maxx, sec.sub_mod, sec.rat_mod, branch_color);
			branch_color += dawg::residue::branch_inc;
			// if gap_overlap is false, clear the upstream buffer of all sequences
			if(!sec.gap_overlap) {
				foreach(details::sequence_data &v, seqs) {
					v.indels.clear();
				}
			}
		}
		foreach(const section &sec, seg) {
			for(wood::data_type::size_type u=1;u<sec.usertree.data().size();++u) {
				const wood::node &n = sec.usertree.data()[u];
				const sequence &ranc = seqs[sec.metatree[u-n.anc]].seq;
				sequence &seq = seqs[sec.metatree[u]].seq;
				sec.evolve(seq, seqs[sec.metatree[u]].indels, n.length,
					branch_color, ranc.begin(), ranc.end(), maxx);
				branch_color += dawg::residue::branch_inc;
			}
			
		}
		// Align segment
		align(aln, seqs, seg.rex);
	}
}

void dawg::details::matic_section::evolve_upstream(
		sequence &child, indel_data &indels, double T, residue::data_type branch_color,
		mutt &m) const {
	double dM, d;
	//insertion and deletion rates
	double ins_rate = ins_mod.rate();
	//double indel_rate = ins_rate+del_rate;
	
	indels.clear();
	
	indel_data::stack del_up, ins_up;

	// Calculate Upstream Deletions
	dM = del_mod.upstream_rate();
	if(dM > DBL_EPSILON) {
		d = m.rand_exp(dM);
		while(d < T) {
			del_up.push(indel_data::element(d/T,
				del_mod.sample_upstream_overlap(m)));
			d += m.rand_exp(dM);
		}
	}
	// Calculate Immortal Link Insertion
	dM = ins_rate;
	if(dM > DBL_EPSILON) {
		d = m.rand_exp(dM);
		while(d < T) {
			ins_up.push(indel_data::element(d/T, ins_mod(m)));
			d += m.rand_exp(dM);
		}
	}
	// Process Immortal Link Insertions
	sequence temp;
	while(!ins_up.empty()) {
		indels.ins.push(ins_up.top());
		ins_up.pop();
		indel_data::element &n = indels.ins.top();
		// Add any upstream indels that happen after
		while(!del_up.empty() && del_up.top().first >= n.first) {
			indels.del.push(del_up.top());
			del_up.pop();
		}
		evolve_indels(child, indels, T, branch_color, temp.begin(),temp.end(),m);
	}

	// Add any outstanding upstream deletions
	while(!del_up.empty()) {
		indels.del.push(del_up.top());
		del_up.pop();
	}
}

dawg::sequence::const_iterator
dawg::details::matic_section::evolve_indels(
		sequence &child, indel_data &indels, double T, residue::data_type branch_color,
		sequence::const_iterator first, sequence::const_iterator last,
		mutt &m) const {
	double f, t;
	double ins_rate = ins_mod.rate(), del_rate = del_mod.rate();
	//double indel_rate = ins_rate+del_rate;

	for(;;) {
		// Is there a deletion that needs to be processed?
		if(!indels.del.empty()) {
			indel_data::element &r = indels.del.top();
			// Something has to be deleted
			if(!indels.ins.empty()) {
				// Deleted Insertion
				indel_data::element &n = indels.ins.top();
				// Did anything happen between the deletion and this insertion
				assert(r.first >= n.first && 1.0 > r.first-n.first );
				t = r.first-n.first;
				double tt = n.first;
				boost::uint32_t u = std::min(r.second, n.second);
				// Determine where the next event occurs
				boost::uint32_t x = next_indel(m.rand_exp(t*T), f, true);
				if(x <= 2*u) {
					// the next event occured between these two
					// how may sites are deleted
					u = x/2;
					// adjust stacks
					r.second -= u;
					n.second -= u;
					if(r.second == 0)
						indels.del.pop();
					if(n.second == 0)
						indels.ins.pop();					
					// push the new event on the proper stack
					if((x&1) == 1) {
						indels.del.push(indel_data::element(tt+t*f/del_rate, del_mod(m)));
					} else {
						do {
							indels.ins.push(indel_data::element(tt+t*f/ins_rate, ins_mod(m)));
							f += m.rand_exp(t*T);
						} while( f < ins_rate );
					}					
				} else {
					r.second -= u;
					n.second -= u;
					if(r.second == 0)
						indels.del.pop();
					if(n.second == 0)
						indels.ins.pop();					
				}
				
				// insert u "deleted insertions" into buffer
				child.insert(child.end(), u, residue(gap_base,
					residue::rate_type(0.0), branch_color));
				// remove u sites from both stacks, pop if empty
			} else if(first != last) {
				// Deleted Original
				assert(r.first < 1.0);
				t = r.first;
				// Determine where the next event occurs
				boost::uint32_t x = next_indel(m.rand_exp(t*T), f, true);
				if(x <= 2*r.second) {
					// copy and mark at most x/2 nucleotides as deleted
					boost::uint32_t u = mark_del(x/2, child, first, last);
					if(u == r.second)
						indels.del.pop();
					else
						r.second -= u;
					if((x&1) == 1) {
						indels.del.push(indel_data::element(t*f/del_rate, del_mod(m)));
					} else {
						do {
							indels.ins.push(indel_data::element(t*f/ins_rate, ins_mod(m)));
							f += m.rand_exp(t*T);
						} while( f < ins_rate );
					}										
				} else {
					// copy and mark at most r.second nucleotides as deleted
					boost::uint32_t u = mark_del(r.second, child, first, last);
					// remove sites from stack, pop if empty
					if(u == r.second)
						indels.del.pop();
					else
						r.second -= u;
				}
			} else {
				// everything possible has been deleted
				break;
			}
		} else if(!indels.ins.empty()) {
			indel_data::element &n = indels.ins.top();
			assert(n.first < 1.0);
			t = 1.0-n.first;
			// Find location of next event
			boost::uint32_t x = next_indel(m.rand_exp(t*T), f, false);
			boost::uint32_t u = n.second;
			if(x <= 2*u) {
				// next event overlaps this one
				// remove u sites from the location and pop if empty
				u = x/2;
				n.second -= u;
				double tt = n.first;
				if(n.second == 0)
					indels.ins.pop();
				// push the new event(s) on the proper stack
				if((x&1) == 1) {
					indels.del.push(indel_data::element(tt+t*f/del_rate, del_mod(m)));
				} else {
					do {
						indels.ins.push(indel_data::element(tt+t*f/ins_rate, ins_mod(m)));
						f += m.rand_exp(t*T);
					} while( f < ins_rate );
				}
			} else {
				indels.ins.pop();
			}
			// Insert u random nucleotides into buffer
			while(u--)
				child.push_back(residue(sub_mod(m),
					static_cast<residue::rate_type>(rat_mod(m)), branch_color));
		} else {
			break; // nothing to do
		}
	}
	return first;
}

void dawg::details::matic_section::evolve(
		sequence &child, indel_data &indels, double T, residue::data_type branch_color,
		sequence::const_iterator first, sequence::const_iterator last,
		mutt &m) const {
		
	//process any existing indels.
	first = evolve_indels(child, indels, T, branch_color, first, last, m);
	
	const double ins_rate = ins_mod.rate(), del_rate = del_mod.rate();
	const double indel_rate = ins_rate+del_rate;
	const double uni_scale = sub_mod.uniform_scale();
	double d = m.rand_exp(T);
	for(;;) {
		sequence::const_iterator start = first;
		// TODO: Optimize out this if?
		// TODO: Variant for constant rate_scale
		// TODO: Optimzie out uni_scale multiplication by changing T and indel_rate
		// TODO: Move to residue_model class
		for(;first != last; ++first) {
			if(first->base() == gap_base)
				continue;
			if(d < indel_rate+first->rate_scalar()*uni_scale)
				break;
			d -= indel_rate+first->rate_scalar()*uni_scale;
		}
		// copy unmodified sites into buffer.
		child.insert(child.end(), start, first);
		if(first == last)
			break;
		if(d < del_rate) {
			indels.del.push(indel_data::element(d/del_rate, del_mod(m)));
			first = evolve_indels(child, indels, T, branch_color, first, last, m);
			d = m.rand_exp(T);
			continue;
		} else
			d -= del_rate;
		double w = first->rate_scalar()*uni_scale;
		residue rez = *first;
		++first;
		while(d < w) {
			rez.base(sub_mod(m,rez.base()));
			d += m.rand_exp(T);
		}
		d -= w;
		child.push_back(rez);
		if(d < ins_rate) {
			do {
				indels.ins.push(indel_data::element(d/ins_rate, ins_mod(m)));
				d += m.rand_exp(T);
			} while(d < ins_rate);
			first = evolve_indels(child, indels, T, branch_color, first, last, m);
		}
		d -= ins_rate;
	}
}

struct aligner_data {
	aligner_data(const sequence &xseq, std::string &xstr) :
		it(xseq.begin()), last(xseq.end()), str(&xstr) {
	}
	sequence::const_iterator it, last;
	std::string *str;
};

void dawg::matic::align(alignment& aln, const seq_buffers_type &seqs, const residue_exchange &rex) {
	assert(aln.size() <= seqs.size());
	
	//unsigned uFlags = 0; //temporary
	
	// construct a table to hold alignment information
	// TODO: Cache this?
	std::vector<aligner_data> aln_table;
	for(alignment::size_type u=0;u<aln.size();++u) {
		aln_table.push_back(aligner_data(seqs[u].seq, aln[u].seq));
	}
	
	// Alignment rules:
	// Insertion & Deleted Insertion  : w/ ins, deleted ins, or gap
	// Deletion & Original Nucleotide : w/ del, original nucl

	unsigned int uStateQuit = rex.is_keep_empty() ? (residue::base_mask+1)*2
		                                          : (residue::base_mask+1)*2-1 ;
	unsigned int uBranch = 0;
	unsigned int uBranchN = 0;
	// Go through each column, adding gaps where neccessary
	for(;;) {
		unsigned int uState =  uStateQuit; // Set to quit
		uBranch = 0; // Set to lowest branch
		// Find column state(s)
		foreach(aligner_data &v, aln_table) {
			if(v.it == v.last)
				continue; // Sequence is done
			uBranchN = v.it->branch();
			if(uBranchN == uBranch) {
				uState &= v.it->base();
			} else if(uBranchN > uBranch) {
				uBranch = uBranchN;
				uState = v.it->base() & uStateQuit;
			}
		}
		switch(((uState+1) >> residue::base_bit_width)&3) {
			case 2: goto ENDFOR; // Yes, you shouldn't use goto, except here
			case 1: // Empty column that we want to ignore
				foreach(aligner_data &v, aln_table) {
					if(v.it != v.last && v.it->branch() == uBranch)
						++v.it;
				}
				break;
			case 0: // Unempty column
				foreach(aligner_data &v, aln_table) {
					if(v.it == v.last || v.it->branch() != uBranch) {
						rex.append_ins(*v.str);
					} else {
						rex.append_residue(*v.str, *(v.it++));
					}
				}
		};
	}
ENDFOR:
	/*noop*/;
}
