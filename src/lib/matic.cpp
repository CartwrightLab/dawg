/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/utils/foreach.h>
#include <dawg/matic.h>
#include <dawg/log.h>
#include <dawg/wood.h>

using namespace dawg;

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
	
	// parse tree and find all named nodes
	info->usertree.parse(ma.tree_tree.begin(), ma.tree_tree.end());
	foreach(const dawg::wood::node &n, info->usertree.data) {
		if(n.label.empty())
			continue;
		if(!info->node_names.insert(n.label).second)
			return DAWG_ERROR("invalid tree; node label '" << n.label
			                  << "' used more than once by Tree.Tree.");
	}
	
	// test whether descendents already exist in this segment
	foreach(section_info &r, seg) {
		std::set<std::string>::const_iterator it = has_intersection(
			r.node_names.begin(), r.node_names.end(),
			info->node_names.begin(), info->node_names.end());
		if(it != info->node_names.end())
			return DAWG_ERROR("invalid tree; descendent '" <<  *it
			               << "' already exists in segment #"
			               << ma.root_segment);
	}
	
	
	// find location to insert
	segment_info::iterator it;
	for(it = seg.begin(); it != seg.end()
	    && info->node_names.count(it->root_name); ++it)
		/*noop*/;
	seg.insert(it, info);
	
	return true;
}

