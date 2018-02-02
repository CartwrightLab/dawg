#pragma once
#ifndef DAWG_MA_H
#define DAWG_MA_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <string>
#include <vector>
#include <iostream>

#include <dawg/trick.h>
#include <dawg/utils/vecio.h>

#include <boost/algorithm/string/case_conv.hpp>

namespace dawg {

///////////////////////////////////////////////////////////
/// \brief dawg::ma is a "model argument" structure
///	std::string subst_model 
///	std::vector<double> subst_params
///	std::vector<double> subst_freqs
///	std::string subst_rate_model
///	std::vector<double> subst_rate_params 

///	std::vector<std::string> indel_model_ins 
///	std::vector<double> indel_params_ins 
///	std::vector<double> indel_rate_ins 
///	unsigned int indel_max_ins 
///	std::vector<std::string> indel_model_del 
///	std::vector<double> indel_params_del 
///	std::vector<double> indel_rate_del 
///	unsigned int indel_max_del 
///
///	std::string tree_model 
///	std::vector<double> tree_params 
///	std::string tree_tree 
///	double tree_scale
///
///	unsigned int root_length 
///	std::string root_seq 
///	std::vector<double> root_rates 
///	unsigned int root_code 
///	unsigned int root_segment 
///	bool root_gapoverlap 
///
///	bool output_markins 
///	bool output_keepempty 
///	bool output_lowercase 
///	bool output_rna 
///////////////////////////////////////////////////////////
struct ma {
#	define XM(name, type, def, desc) type XV(name) ;
#	include <dawg/details/dawgma.xmh>
#	undef XM
	std::string name;

	ma(const std::string &_n = std::string() ) :
#	define XM(name, type, def, desc) XV(name) (def),
#	include <dawg/details/dawgma.xmh>
#	undef XM
	name(_n)
	{ }
	
	static bool from_trick(const dawg::trick &trk, std::vector<dawg::ma> &v);
	
	void read_section(const trick::data_type::value_type &sec);
	
	template<class CharType, class CharTrait>
	static void help(std::basic_ostream<CharType, CharTrait>& o);
	
private:
};

template<class CharType, class CharTrait>
inline std::basic_ostream<CharType, CharTrait>& 
operator<<(std::basic_ostream<CharType, CharTrait>& o, const ma &a) {
	if(!o.good()) return o;
		
	o << set_open('\x7f') << set_close('\x7f') << set_delimiter(',');

	o << "[[ " << a.name << " ]]" << std::endl;
#	define XM(name, type, def, desc) o << XP(name) " = " << a.XV(name) << std::endl;
#	include <dawg/details/dawgma.xmh>
#	undef XM

	return o;
}

namespace details {
inline std::string ma_help_name(const char *cs) {
	std::string ret(cs);
	typedef boost::iterator_range<std::string::iterator> range;
	range r(ret.begin(), ret.begin()+1);
	boost::to_upper(r);
	std::size_t pos = 0;
	while((pos = ret.find('.', pos+1)) != std::string::npos) {
		r = range(ret.begin()+pos+1, ret.begin()+pos+2);
		boost::to_upper(r);
	}
	return ret;
}
};

template<class CharType, class CharTrait>
void ma::help(std::basic_ostream<CharType, CharTrait>& o) {
	o << "[REGULAR PARAMETERS]" << std::endl;
#	define XM(name, type, def, desc) o << details::ma_help_name(XP(name)) \
	<< " - " << (desc) << std::endl;
#	include <dawg/details/dawgma.xmh>
#	undef XM
	o << "\n[GLOBAL PARAMETERS]" << std::endl;
#	define XM(name, type, def, desc) o << details::ma_help_name(XP(name)) \
	<< " - " << (desc) << std::endl;
#	include <dawg/details/global.xmh>
#	undef XM
	
}


} //namespace dawg
#endif //DAWG_MA_H
