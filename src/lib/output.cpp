/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/log.h>

#include <dawg/output.h>
#include <dawg/utils/foreach.h>

#include <cstring>
#include <iostream>

#include <boost/range/iterator_range.hpp>

using namespace dawg;
using namespace std;

bool dawg::output::open(const char *file_name) {
	typedef boost::iterator_range<const char*> cs_range;

	format_id = 0;
	cs_range format(file_name, file_name);
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
		} else if((mid = strrchr(file_name, '.')) != NULL) {
			// file.format
			format = boost::make_iterator_range(mid+1, strchr(mid+1, '\0'));
		}
		if(format && !set_format(format)) {
			return DAWG_ERROR("unknown output format \'"
						<< std::string(format.begin(), format.end()) << "\'.");
		}
	}
	if(file_name != NULL && file_name[0] != '\0'
		&& strcmp(file_name, "-") != 0) {
		if(fout.is_open())
			fout.close();
		fout.open(file_name, ios_base::out|ios_base::trunc);
		if(!fout.is_open()) {
			return DAWG_ERROR("unable to open output file \'" << file_name<< "\'.");
		}
		set_ostream(fout);
	} else {
		set_ostream(cout);
	}
	
	do_op = &output::print_fasta;

	return true;
}

void dawg::output::print_fasta(const alignment& aln) {
	ostream &out = *p_out;
	
	foreach(const alignment::value_type& v, aln) {
		out << ">" << v.label << endl;

		const char *it = v.seq.c_str();
		const char *last = it+v.seq.size();
		for(;it+75 < last;it+=75) {
			out.write(it, 75);
			out << endl;
		}
		out.write(it, last-it);
		out << endl << endl;
	}
}
