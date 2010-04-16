/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/log.h>

#include <dawg/output.h>
#include <dawg/utils.h>

#include <cstring>
#include <utility>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/range/iterator_range.hpp>

using namespace dawg;
using namespace std;

bool dawg::output::open(const char *file_name) {
	typedef boost::iterator_range<const char*> cs_range;

	static const char format_keys[][10] = {
		"aln", "poo", "fasta", "fas", "fsa",
		"nexus", "nex", "phylip", "phy"
	};
	std::size_t format_id;
	cs_range format;
	if(file_name != NULL && file_name[0] != '\0') {
		const char *mid = strchr(file_name, ':');
#ifdef BOOST_WINDOWS
		if(mid != NULL && mid != file_name+1) {
#else
		if(mid != NULL) {
#endif		
			// format:file
			format = boost::make_iterator_range(file_name, mid);
			file_name = mid+1;
		} else {
			// file.format
			mid = strrchr(file_name, '.');
			format = boost::make_iterator_range(mid+1, strchr(mid+1, '\0'));
		}
		if(!format) {
			format_id = key_switch(format, format_keys);
			if(format_id == -1) {
				return DAWG_ERROR("unknown output format \'"
				       << std::string(format.begin(), format.end()) << "\'.");
			}
		}
	}
	if(file_name != NULL && file_name[0] != '\0'
		&& strcmp(file_name, "-") != 0) {
		fout.close();
		fout.open(file_name, ios_base::out|ios_base::trunc);
		if(!fout.is_open()) {
			return DAWG_ERROR("unable to open output file \'" << file_name<< "\'.");
		}
	}
	p_out = fout.is_open() ? static_cast<ostream*>(&fout)
	                       : static_cast<ostream*>(&cout);
	                       
	*p_out << "Hello World" << endl
	       << format_id << " " << std::string(format.begin(), format.end()) << endl;
	return true;
}


