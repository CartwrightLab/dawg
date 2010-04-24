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
#include <boost/spirit/include/karma_generate.hpp>
#include <boost/spirit/include/karma_uint.hpp>
#include <boost/spirit/include/karma_right_alignment.hpp>
#include <boost/spirit/include/karma_char_.hpp>

using namespace dawg;
using namespace std;

bool dawg::output::open(const char *file_name, unsigned int max_rep, bool append) {
	typedef boost::iterator_range<const char*> cs_range;
	set_format("aln");
	rep = 0;
	cs_range format(file_name, file_name);
	const char *mid = NULL;
	if(file_name != NULL && file_name[0] != '\0') {
		mid = strchr(file_name, ':');
#ifdef BOOST_WINDOWS
		if(mid != NULL && mid != file_name+1) {
#else
		if(mid != NULL) {
#endif		
			// format:file
			format = boost::make_iterator_range(file_name, mid);
			file_name = mid+1;
			// find extension point for later
			mid = strrchr(file_name, '.');
		} else if((mid = strrchr(file_name, '.')) != NULL) {
			// file.format
			format = boost::make_iterator_range(mid+1, strchr(mid+1, '\0'));
		}
		if(format && !set_format(format)) {
			return DAWG_ERROR("unknown output format \'"
						<< std::string(format.begin(), format.end()) << "\'.");
		}
	}
	if(file_name != NULL && file_name[0] != '\0' && strcmp(file_name, "-") != 0) {
		// set append option before we open the file
		app = append;
		// find how wide the id number must be
		unsigned int m = max_rep;
		split_width = (m == 0) ? 0 : (m < 10) ? 1 : (m < 100) ? 2 : (m < 1000) ? 3 :
			(m < 10000) ? 4 : (m < 100000) ? 5 : (m < 1000000) ? 6 :
			(m < 10000000) ? 7 : (m < 100000000) ? 8 :
			(m < 1000000000) ? 9 : 10;
		// open our omnibus output if desired
		if(split_width == 0) {
			if(!open_file(file_name))
				return DAWG_ERROR("unable to open output file \'" << file_name<< "\'.");
		} else {
			// setup output_filename
			if(mid == NULL) {
				split_file_name.assign(file_name);
				split_file_name.append(1, '-');
				split_file_name.append(split_width, '0');
				split_id_offset = split_file_name.size()-split_width;
			} else {
				split_file_name.assign(file_name, mid);
				split_file_name.append(1, '-');
				split_file_name.append(split_width, '0');
				split_file_name.append(mid);
				split_id_offset = (mid-file_name)+1;
			}
		}
	} else {
		// turn off appending and spliting
		app = false;
		split_width = 0;
		set_ostream(cout);
	}
	return true;
}

bool dawg::output::open_file(const char* file_name) {
	if(fout.is_open())
		fout.close();
	ios_base::openmode om = ios_base::out | (app ? ios_base::app : ios_base::trunc);
	fout.open(file_name, om);
	if(!fout.is_open()) {
		set_ostream(NULL);
		return false;
	}
	set_ostream(fout);
	return true;
}

bool dawg::output::open_next() {
	// no need to open next file
	if(split_width == 0)
		return true;
	using boost::spirit::karma::generate;
	using boost::spirit::karma::uint_;
	using boost::spirit::karma::lit;
	using boost::spirit::karma::right_align;

	string::iterator it = split_file_name.begin()+split_id_offset;
	generate(it, right_align(split_width,lit('0'))[uint_], rep);
	
	if(!open_file(split_file_name.c_str()))
		return DAWG_ERROR("unable to open output file \'" << split_file_name << "\'.");
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
		<< PACKAGE_STRING << ")\n" << endl;

	// find alignment length
	string::size_type u=0, len = aln[0].seq.length();
	for(;u+60 < len;u += 60) {
		foreach(const alignment::value_type& v, aln) {
			out << setw(aln.max_label_width_14) << left << v.label << ' ';
			out.write(&v.seq[u], 60);
			out << endl;
		}
		out << endl;
	}
	len = len-u;
	foreach(const alignment::value_type& v, aln) {
		out << setw(aln.max_label_width_14) << left << v.label << ' ';
		out.write(&v.seq[u], len);
		out << endl;
	}
	out << endl;
}

void dawg::output::print_poo(const alignment& aln) {
	ostream &out = *p_out;

	foreach(const alignment::value_type& v, aln) {
		out << setw(aln.max_label_width) << v.label
			<< ' ' << v.seq
			<< endl;
	}
	out << endl;
}

void dawg::output::print_phylip(const alignment& aln) {
	ostream &out = *p_out;

	out << "  " << aln.size() << "    " << aln[0].seq.size() << endl;
	foreach(const alignment::value_type& v, aln) {
		out << setw(10) << left << v.label.substr(0,10) << v.seq << endl;
	}
	out << endl;
}

// TODO: save sequence type in aln
void dawg::output::print_nexus(const alignment& aln) {
	ostream &out = *p_out;
	
	static char datatypes[][10] = {
		"DNA", "RNA", "PROTEIN", "DNA",
		"DNA", "RNA", "PROTEIN", "DNA"
	};
	
	// TODO: Move this to the block support routine
	out << "#NEXUS\n[Created by " PACKAGE_STRING "]" << endl;

	out << "BEGIN DATA;\n"
		   "\tDIMENSIONS NTAX=" << aln.size() << " NCHAR=" << aln[0].seq.size() << ";\n"
	       "\tFORMAT DATATYPE=";
	
	out << datatypes[aln.seq_type];
	out << " MISSING=? GAP=- MATCHCHAR=. EQUATE=\"+=-\";\n"
	       "\tMATRIX" << endl;
	
	// Write sequences in non-interleaved format
	foreach(const alignment::value_type& v, aln) {
		out << setw(aln.max_label_width_14) << left << v.label << ' ' << v.seq << endl;
	}
	// Close data block
	out << ";\nEND;\n" << endl;

}
