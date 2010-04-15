/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/residue.h>
#include <dawg/matic.h>
#include <dawg/log.h>
#include <dawg/wood.h>

#include <dawg/utils/foreach.h>
#include <dawg/utils/vecio.h>

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
	if(!info->sub_mod.create(ma.subst_model,
		ma.subst_params.begin(), ma.subst_params.end(),
		ma.subst_freqs.begin(), ma.subst_freqs.end()))
		return DAWG_ERROR("substitution model could not be created.");
	
	if(seg.empty()) { // new segment
		seg.rex.model(info->sub_mod.seq_type(), ma.output_lowercase,
			ma.output_markins, ma.output_keepempty);
	} else if(!seg.rex.is_same_model(info->sub_mod.seq_type(),
	           ma.output_lowercase, ma.output_markins, ma.output_keepempty)) {
		return DAWG_ERROR("the sequence type or format options of a section is different than its segment.");
	}	
	
	if(!info->rat_mod.create(ma.subst_rate_model, ma.subst_rate_params.begin(),
		ma.subst_rate_params.end()))
		return DAWG_ERROR("heterogenous rate model could not be created.");
	if(!info->ins_mod.create(ma.indel_model_ins.begin(), ma.indel_model_ins.end(),
		ma.indel_rate_ins.begin(), ma.indel_rate_ins.end(),
		ma.indel_params_ins.begin(), ma.indel_params_ins.end()))
		return DAWG_ERROR("insertion model could not be created.");
	if(!info->del_mod.create(ma.indel_model_del.begin(), ma.indel_model_del.end(),
		ma.indel_rate_del.begin(), ma.indel_rate_del.end(),
		ma.indel_params_del.begin(), ma.indel_params_del.end()))
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
	foreach(label_to_index_type::value_type &lab, label_union) {
		lab.second = u++;
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


void dawg::matic::walk(alignment& aln) {
	if(configs.empty())
		return;
	aln.resize(label_union.size());
	vector<details::sequence_data> seqs(label_union.size());
	int uu = 0;
	foreach(const label_to_index_type::value_type &kv, label_union) {
		aln[uu].seq.clear();
		seqs[uu].indels.clear();
		seqs[uu].seq.clear();
		//TODO: Move this to pre walk?
		aln[uu++].label = kv.first;
	}
	
	{	// take first seg and do upstream indels
		// root.gapoverlap = false prevents upstream indel creation
		const segment &seg = configs[0];
		branch_color = 0;
		foreach(const section &sec, seg) {
			if(!sec.gap_overlap)
				continue;
			for(wood::data_type::size_type u=1;u<sec.usertree.data().size();++u) {
				const wood::node &n = sec.usertree.data()[u];
				details::sequence_data &sd = seqs[sec.metatree[u]];
				sec.evolve_upstream(sd.seq, sd.indels, n.length,
					branch_color, maxx);
			}
		}
		align(aln, seqs, seg.rex);
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
		}		
		foreach(const section &sec, seg) {
			if(!sec.gap_overlap) {
				// if gap_overlap is false, clear the upstream buffer
				// of all sequences in the tree
				for(wood::data_type::size_type u=1;u<sec.usertree.data().size();++u) {
					seqs[sec.metatree[u]].indels.clear();
				}
			}
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
	double ins_rate = ins_mod.rate(), del_rate = del_mod.rate();
	double indel_rate = ins_rate+del_rate;
	
	indels.clear();
	
	indel_data::stack del_up;

	// Calculate Upstream Deletions
	dM = del_rate*(del_mod.meansize()-1.0);
	if(dM > DBL_EPSILON) {
		d = m.rand_exp(dM);
		while(d < T) {
			boost::uint32_t u;
			do {
				u = del_mod(m);
			} while(u == 1);
			del_up.push(indel_data::element(d/T,  1+m.rand_uint32(u-1)));
			d +=m.rand_exp(dM);
		}
	}
	// Calculate Immortal Link Insertions
	if(ins_rate > DBL_EPSILON) {
		dM = ins_rate;
		d = m.rand_exp(dM);
		while(d < T) {
			indels.ins.push(indel_data::element(d/T, ins_mod(m)));
			d += m.rand_exp(dM);
		}
	}
	// Process any Imortal Link Insertions
	while(!indels.ins.empty()) {
		// Fetch most recent immortal link insertion
		indel_data::element &n = indels.ins.top();
		// Did any upstream deletions occur before this?
		while(!del_up.empty() && del_up.top().first >= n.first) {
			indels.del.push(del_up.top());
			del_up.pop();
		}
		evolve_indels(child, indels, T, branch_color,
			sequence::const_iterator(), sequence::const_iterator(), m);
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
	double indel_rate = ins_rate+del_rate;

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
				boost::uint32_t u = min(r.second, n.second);
				// Determine where the next event occurs
				boost::uint32_t x = next_indel(m.rand_exp(t*T), f);
				if(x < 2*u) {
					// the next event occured between these two
					// how may sites are deleted
					u = (x+1)/2;
					// push the new event on the proper stack
					if((x&1) == 1) {
						indels.ins.push(indel_data::element(n.first+t*f/ins_rate, ins_mod(m)));
					} else {
						indels.del.push(indel_data::element(n.first+t*f/del_rate, del_mod(m)));
					}
				}
				// insert u "deleted insertions" into buffer
				child.insert(child.end(), u, residue(residue::deleted_base,
					residue::rate_type(1.0), branch_color));
				// remove u sites from both stacks, pop if empty
				r.second -= u;
				n.second -= u;
				if(r.second == 0)
					indels.del.pop();
				if(n.second == 0)
					indels.ins.pop();
			} else if(first != last) {
				// Deleted Original
				assert(r.first < 1.0);
				t = r.first;
				boost::uint32_t u = r.second;
				// Determine where the next event occurs
				boost::uint32_t x = next_indel(m.rand_exp(t*T), f);
				if(x < 2*u) {
					// the next event occured between these two
					// how may sites are deleted
					u = (x+1)/2;
					// push the new event on the proper stack
					if((x&1) == 1) {
						indels.ins.push(indel_data::element(t*f/ins_rate, ins_mod(m)));
					} else {
						indels.del.push(indel_data::element(t*f/del_rate, del_mod(m)));
					}
				}
				// copy and mark at most u nucleotides as deleted
				boost::uint32_t uu;
				for(uu=0;uu != u && first != last;++first) {
					child.push_back(*first);
					if(first->is_deleted())
						continue;
					child.back().base(residue::deleted_base);
					++uu;
				}
				// remove sites from stack, pop if empty
				r.second -= uu;
				if(r.second == 0)
					indels.del.pop();
			} else {
				// everything possible has been deleted
				break;
			}
		} else if(!indels.ins.empty()) {
			indel_data::element &n = indels.ins.top();
			assert(n.first < 1.0);
			t = 1.0-n.first;
			// Find location of next event
			boost::uint32_t x = next_indel(m.rand_exp(t*T), f)-1;
			boost::uint32_t u = n.second;
			if(x <= 2*u) {
				// next event overlaps this one
				u = x/2;
				// push the new event on the proper stack
				if((x&1) == 0) {
					indels.ins.push(indel_data::element(n.first+t*f/ins_rate, ins_mod(m)));
				} else {
					indels.del.push(indel_data::element(n.first+t*f/del_rate, del_mod(m)));
				}
			}
			// remove u sites from the location and pop if empty
			n.second -= u;
			if(n.second == 0)
				indels.ins.pop();
			// Insert u random nucleotides into buffer
			while(u--)
				child.push_back(residue(sub_mod(m),
					static_cast<residue::rate_type>(rat_mod(m)), branch_color));
		} else {
			// nothing to do
			break;
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
	
	double ins_rate = ins_mod.rate(), del_rate = del_mod.rate();
	double indel_rate = ins_rate+del_rate;	
	double d = m.rand_exp(T);
	for(;;) {
		sequence::const_iterator start = first;
		// TODO: Optimize out this if?
		// TODO: Variant for constant rate_scale
		for(;first != last && first->rate_scalar()+indel_rate <= d; ++first) {
			if(!first->is_deleted())
				d -= indel_rate+first->rate_scalar();
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
		double w = first->rate_scalar();
		residue rez = *first;
		++first;
		while(d < w) {
			rez.base(sub_mod(m,rez.base()));
			// how much space is left in the substitution section
			w = w - d;
			d = m.rand_exp(T);		
		}
		d -= w;
		child.push_back(rez);
		if(d < ins_rate) {
			indels.ins.push(indel_data::element(d/ins_rate, ins_mod(m)));
			first = evolve_indels(child, indels, T, branch_color, first, last, m);
			d = m.rand_exp(T);
		} else
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
	
	unsigned uFlags = 0; //temporary
	
	// construct a table to hold alignment information
	// TODO: Cache this?
	std::vector<aligner_data> aln_table;
	for(alignment::size_type u=0;u<aln.size();++u) {
		aln_table.push_back(aligner_data(seqs[u].seq, aln[u].seq));
	}
	
	// Alignment rules:
	// Insertion & Deleted Insertion  : w/ ins, deleted ins, or gap
	// Deletion & Original Nucleotide : w/ del, original nucl

	unsigned int uStateQuit = rex.is_keep_empty() ? (residue::deleted_base+1)*2
		                                          : (residue::deleted_base+1)*2-1 ;
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
						residue_exchange::str_type s = rex.decode_ins();
						v.str->append(s.first, s.second);
					} else {
						residue_exchange::str_type s = rex.decode(*v.it);
						v.str->append(s.first, s.second);
						++v.it;
					}
				}
		};
	}
ENDFOR:
	/*noop*/;
}

