/****************************************************************************
 *  Copyright (C) 2010-2012 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/
#include <boost/range/iterator_range.hpp>
#include <boost/spirit/include/karma_char_.hpp>
#include <boost/spirit/include/karma_generate.hpp>
#include <boost/spirit/include/karma_uint.hpp>
#include <boost/spirit/include/karma_right_alignment.hpp>

#include <cstring>
#include <iostream>
#include <iomanip>

#include <dawg/details/config.h>
#include <dawg/log.h>
#include <dawg/error.h>

#include <dawg/output.h>

using namespace dawg;
using namespace std;

bool dawg::output::open(const char *file_name, unsigned int max_rep,
	bool split, bool append, bool label) {
	typedef boost::iterator_range<const char*> cs_range;
	set_format("aln");
	last_rep = max_rep;
	rep = 0;
	cs_range format(file_name, file_name);
	const char *mid = nullptr;
	if(file_name != nullptr && file_name[0] != '\0') {
		mid = strchr(file_name, ':');
#ifdef BOOST_WINDOWS
		if(mid != nullptr && mid != file_name+1) {
#else
		if(mid != nullptr) {
#endif		
			// format:file
			format = boost::make_iterator_range(file_name, mid);
			file_name = mid+1;
			// find extension point for later
			mid = strrchr(file_name, '.');
		} else if((mid = strrchr(file_name, '.')) != nullptr) {
			// file.format
			format = boost::make_iterator_range(mid+1, (const char*)strchr(mid+1, '\0'));
		}
		if(format && !set_format(format)) {
			throw dawg::dawg_error_t(dawg_error::unknown_output_format, \
			    std::string(format.begin(), format.end()));
		}
	}
	label_width = 1+static_cast<unsigned int>(log10(1.0*max_rep));
	current_label.assign(label_width, '0');
	do_label = label;
	
	if(file_name != nullptr && file_name[0] != '\0' && strcmp(file_name, "-") != 0) {
		// set append and split options before we open the file
		do_append = append;
		do_split  = split;
		// open our omnibus output if desired
		if(!do_split) {
			if(!open_file(file_name)) {
				throw dawg::dawg_error_t(dawg_error::open_output_file_fail, \
				    std::string(file_name));
			}
		} else {
			// setup output_filename
			if(mid == nullptr) {
				split_file_name.assign(file_name);
				split_file_name.append(1, '-');
				split_file_name.append(current_label);
				split_id_offset = split_file_name.size()-label_width;
			} else {
				split_file_name.assign(file_name, mid);
				split_file_name.append(1, '-');
				split_file_name.append(current_label);
				split_file_name.append(mid);
				split_id_offset = (mid-file_name)+1;
			}
		}
	} else {
		// turn off appending and spliting
		do_append = false;
		do_split = false;
		set_ostream(cout);
	}
	return true;
}

bool dawg::output::open_file(const char* file_name) {
	if(fout.is_open())
		fout.close();
	ios_base::openmode om = ios_base::out |
		(do_append ? ios_base::app : ios_base::trunc);
	fout.open(file_name, om);
	if(!fout.is_open()) {
		set_ostream(nullptr);
		return false;
	}
	set_ostream(fout);
	return true;
}

bool dawg::output::open_next() {
	// Generate next id number
	using boost::spirit::karma::generate;
	using boost::spirit::karma::uint_;
	using boost::spirit::karma::lit;
	using boost::spirit::karma::right_align;
	generate(current_label.begin(), right_align(label_width, lit('0'))[uint_], rep);

	// no need to open next file
	if(!do_split)
		return true;

	// replace 
	split_file_name.replace(split_id_offset, label_width, current_label);
	
	if(!open_file(split_file_name.c_str())) {
		throw dawg::dawg_error_t(dawg_error::open_input_file_fail, std::string(split_file_name));
	}
	return true;
}

void dawg::output::print_poo(const alignment& aln) {
	ostream &out = *p_out;

	for(const alignment::value_type& v : aln) {
		out << setw(aln.max_label_width) << v.label
			<< '_' << current_label
			<< ' ' << v.seq
			<< endl;
	}
	out << endl;
}

void dawg::output::print_fasta(const alignment& aln) {
	ostream &out = *p_out;

	for(const alignment::value_type& v : aln) {
		out << ">" << v.label;
		if(do_label)
			out << "_" << current_label;
		out << endl;

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
	string::size_type aln_label_width = do_label ?
		std::max(aln.max_label_width+label_width+1,string::size_type(14)) :
		std::max(aln.max_label_width,string::size_type(14)) ;
	for(;u+60 < len;u += 60) {
		for(const alignment::value_type& v : aln) {
			out << v.label;
			if(do_label) {
				out << '_' << current_label;
				out << setw(aln_label_width-(v.label.length()+label_width+1)+1) << ' ';
			} else {
				out << setw(aln_label_width-v.label.length()+1) << ' ';
			}
			out.write(&v.seq[u], 60);
			out << endl;
		}
		out << endl;
	}
	len = len-u;
	for(const alignment::value_type& v : aln) {
		out << v.label;
		if(do_label) {
			out << '_' << current_label;
			out << setw(aln_label_width-(v.label.length()+label_width+1)+1) << ' ';
		} else {
			out << setw(aln_label_width-v.label.length()+1) << ' ';
		}
		out.write(&v.seq[u], len);
		out << endl;
	}
	out << endl;
}

void dawg::output::print_phylip(const alignment& aln) {
	ostream &out = *p_out;

	string::size_type aln_label_width = do_label ?
		std::max(aln.max_label_width+label_width+1,string::size_type(9)) :
		std::max(aln.max_label_width,string::size_type(9)) ;

	out << "  " << aln.size() << "    " << aln[0].seq.size() << endl;
	for(const alignment::value_type& v : aln) {
		out << v.label;
		if(do_label) {
			out << '_' << current_label;
			out << setw(aln_label_width-(v.label.length()+label_width+1)+1) << ' ';
		} else {
			out << setw(aln_label_width-v.label.length()+1) << ' ';
		}
		out << v.seq << endl;
	}
	out << endl;
}

// TODO: save sequence type in aln
void dawg::output::print_nexus(const alignment& aln) {
	ostream &out = *p_out;
	
	static constexpr char datatypes[][10] = {
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
	// string::size_type aln_label_width = do_label ?
	// 	aln.max_label_width+label_width+1 :
	// 	aln.max_label_width;
	
	for(const alignment::value_type& v : aln) {
		out << v.label;
		if(do_label) {
			out << '_' << current_label;
		}
		out << setw(aln.max_label_width-v.label.length()+1) << ' ';
		out << v.seq << endl;
	}
	// Close data block
	out << ";\nEND;\n" << endl;

}
