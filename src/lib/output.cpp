/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <dawg/log.h>

#include <dawg/output.h>
#include <dawg/utils/foreach.h>

#include <cstring>
#include <iostream>
#include <iomanip>

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

void dawg::output::print_aln(const alignment& aln) {
	ostream &out = *p_out;
	out << "CLUSTAL multiple sequence alignment (Created by "
		<< PACKAGE_STRING << ")\n\n" << endl;

	// TODO: Cache width
	string::size_type max_width = 14;
	foreach(const alignment::value_type& v, aln) {
		max_width = std::max(max_width, v.label.length());
	}
	// find alignment length
	string::size_type u=0, len = aln[0].seq.length();
	for(;u+60 < len;u += 60) {
		foreach(const alignment::value_type& v, aln) {
			out << setw(max_width) << left << v.label << ' ';
			out.write(&v.seq[u], 60);
			out << endl;
		}
		out << '\n' << endl;
	}
	len = len-u;
	foreach(const alignment::value_type& v, aln) {
		out << setw(max_width) << left << v.label << ' ';
		out.write(&v.seq[u], len);
		out << endl;
	}
	out << '\n' << endl;
}

void dawg::output::print_poo(const alignment& aln) {
	ostream &out = *p_out;

	// TODO: Cache width
	string::size_type max_width = 0;
	foreach(const alignment::value_type& v, aln) {
		max_width = std::max(max_width, v.label.length());
	}
	foreach(const alignment::value_type& v, aln) {
		out << setw(max_width) << v.label
			<< ' ' << v.seq
			<< '\n' << endl;
	}
	out << endl;
}

void dawg::output::print_phylip(const alignment& aln) {
	ostream &out = *p_out;

	out << "  " << aln.size() << "    " << aln[0].seq.size() << endl;
	foreach(const alignment::value_type& v, aln) {
		out << setw(10) << left << v.label << v.seq << endl;
	}
	out << endl;
}

// TODO: save sequence type in aln
void dawg::output::print_nexus(const alignment& aln) {
	ostream &out = *p_out;

	// TODO: Move this to the block support routine
	out << "#NEXUS\n[Created by " PACKAGE_STRING "]" << endl;

	out << "BEGIN DATA;\n"
		   "\tDIMENSIONS NTAX=" << aln.size() << " NCHAR=" << aln[0].seq.size() << ";\n"
	       "\tFORMAT DATATYPE=";
	out << "DNA"; //TODO or protein
	out << " MISSING=? GAP=- MATCHCHAR=. EQUATE=\"+=-\";\n"
	       "\tMATRIX" << endl;
	
	// Write sequences in non-interleaved format
	foreach(const alignment::value_type& v, aln) {
		out << setw(15) << left << v.label << ' ' << v.seq << endl;
	}
	// Close data block
	out << ";\nEND;\n" << endl;

}



