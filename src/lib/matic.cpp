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
	info->usertree.parse(ma.tree_tree.begin(), ma.tree_tree.end());
		
	// test whether descendents already exist in this segment
	foreach(section &r, seg) {
		std::set<std::string>::const_iterator it = has_intersection(
			r.usertree.descs().begin(), r.usertree.descs().end(),
			info->usertree.descs().begin(), info->usertree.descs().end());
		if(it != info->usertree.descs().end())
			return DAWG_ERROR("invalid tree; descendent '" <<  *it
			               << "' already exists in segment #"
			               << ma.root_segment);
	}
		
	// find location to insert
	segment::iterator it;
	for(it = seg.begin(); it != seg.end()
	    && !info->usertree.has_desc(it->usertree.root_label()); ++it)
		/*noop*/;
	seg.insert(it, info);
	
	return true;
}

typedef map<string,sequence> seq_map;

void dawg::matic::walk() {
	sequence seq_buf, seq_buf2;
	rex.model(residue_exchange::DNA);
	foreach(const segment &seg, configs) {
		seq_map seqs;
		foreach(const section &sec, seg) {
			pair<seq_map::iterator,bool> res =
				seqs.insert(make_pair(sec.usertree.root_label(), sequence()));
			if(res.second)
				sec.rut_mod(res.first->second, maxx, sec.sub_mod, sec.rat_mod);
			wood::data_type::const_iterator nit = sec.usertree.data().begin();
			for(++nit;nit!=sec.usertree.data().end();++nit) {
				seqs[nit->label] = seqs[(nit-nit->anc)->label];
			}
		}
		
		foreach(seq_map::value_type &kv, seqs) {
			string ss(kv.second.size(), ' ');
			rex.decode_array(kv.second.begin(), kv.second.end(), ss.begin());
			cout << kv.first << " "
				 << set_open('\x7f') << set_delimiter('\x7f') << set_close('\x7f')
				 << ss << endl;
		}
		cout << endl;
	}
}

void dawg::details::matic_section::evolve_upstream(
		sequence &child, indel_data &indels, double T,
		sequence::const_iterator first, sequence::const_iterator last) {
	double dM, d;
	//insertion and deletion rates
	double ins_rate = ins_mod.rate(), del_rate = del_mod.rate();
	double indel_rate = ins_rate+del_rate;

	indel_data::stack del_up;

	// Calculate Upstream Deletions
	dM = del_rate*(del_mod.meansize()-1.0);
	if(dM > DBL_EPSILON) {
		d = rand_exp(dM);
		while(d < T) {
			size_type u;
			do {
				u = ();
			} while(u == 1);
			del_up.push(indel_data::element(d/T,  1+rand_uint(u-2)));
			d += rand_exp(dM);
		}
	}
	// Calculate Immortal Link Insertions
	if(ins_rate > DBL_EPSILON) {
		dM = ins_rate;
		d = rand_exp(dM);
		while(d < T) {
			indels.ins.push(indel_data::element(d/T, ins_mod()));
			d += rand_exp(dM);
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
		evolve_indels(seq, indels, dTime, seq.begin(), seq.begin(), sub_mod, ins_mod, del_mod);
	}
	// Add any outstanding upstream deletions
	while(!del_up.empty()) {
		indels.del.push(del_up.top());
		del_up.pop();
	}
}

dawg::sequence::const_iterator
dawg::details::matic_section::evolve_indels(
		sequence &child, indel_data &indels, double T,
		sequence::const_iterator first, sequence::const_iterator last)	double f, t;
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
				t = r.first-n.first;
				boost::uint32_t u = min(r.second, n.second);
				// Determine where the next event occurs
				boost::uint32_t x = next_indel(rand_exp(t*T), f, ins_rate, del_rate);
				if(x < 2*u) {
					// the next event occured between these two
					// how may sites are deleted
					u = (x+1)/2;
					// push the new event on the proper stack
					if((x&1) == 1) {
						indels.ins.push(indel_data::element(n.first+t*f/ins_rate/T, ins_mod()));
					} else {
						indels.del.push(indel_data::element(n.first+t*f/del_rate/T, del_mod()));
					}
				}
				// insert u "deleted insertions" into buffer
				seq.insert(seq.end(), u, residue(0, residue::rate_type(1.0), branch_color, true));
				// remove u sites from both stacks, pop if empty
				r.second -= u;
				n.second -= u;
				if(r.second == 0)
					indels.del.pop();
				if(n.second == 0)
					indels.ins.pop();
			} else if(first != last) {
				// Deleted Original
				t = r.first;
				boost::uint32_t u = r.second;
				// Determine where the next event occurs
				boost::uint32_t x = next_indel(rand_exp(t*T), f, ins_rate, del_rate);
				if(x < 2*u) {
					// the next event occured between these two
					// how may sites are deleted
					u = (x+1)/2;
					// push the new event on the proper stack
					if((x&1) == 1) {
						indels.ins.push(indel_data::element(t*f/ins_rate/T, ins_mod()));
					} else {
						indels.del.push(indel_data::element(t*f/del_rate/T, del_mod()));
					}
				}
				// copy and mark at most u nucleotides as deleted
				boost::uint32_t uu;
				for(uu=0;uu != u && first != last;++first) {
					seq.push_back(*first);
					if(first->is_deleted())
						continue;
					seq.back().mark_deleted(true);
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
			assert(n.first < T);
			t = T-n.first;
			// Find location of next event
			boost::uint32_t x = next_indel(rand_exp(1.0/t), f,ins_rate,del_rate)-1;
			boost::uint32_t u = n.second;
			if(x <= 2*u) {
				// next event overlaps this one
				u = x/2;
				// push the new event on the proper stack
				if((x&1) == 0) {
					indels.ins.push(indel_data::element(n.first+t*f/ins_rate/T, ins_mod()));
				} else {
					indels.del.push(indel_data::element(n.first+t*f/del_rate/T, del_mod()));
				}
			}
			// remove u sites from the location and pop if empty
			n.second -= u;
			if(n.second == 0)
				indels.ins.pop();
			// Insert u random nucleotides into buffer
			while(u--)
				seq.push_back(sub_mod());
		} else {
			// nothing to do
			break;
		}
	}
	return first;
}

dawg::details::matic_section::evolve(
		sequence &child, indel_data &indels, double T,
		sequence::const_iterator first, sequence::const_iterator last) {
}




