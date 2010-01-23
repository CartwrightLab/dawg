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
	segment_info &seg = configs[ma.root_segment];

	// construction section_info
	std::auto_ptr<section_info> info(new section_info);
	
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
	foreach(section_info &r, seg) {
		std::set<std::string>::const_iterator it = has_intersection(
			r.usertree.descs().begin(), r.usertree.descs().end(),
			info->usertree.descs().begin(), info->usertree.descs().end());
		if(it != info->usertree.descs().end())
			return DAWG_ERROR("invalid tree; descendent '" <<  *it
			               << "' already exists in segment #"
			               << ma.root_segment);
	}
		
	// find location to insert
	segment_info::iterator it;
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
	foreach(const segment_info &seg, configs) {
		seq_map seqs;
		foreach(const section_info &sec, seg) {
			pair<seq_map::iterator,bool> res = seqs.insert(make_pair(sec.usertree.root_label(), sequence()));
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

