#pragma once
#ifndef DAWG_OUTPUT_H
#define DAWG_OUTPUT_H
/****************************************************************************
 *  Copyright (C) 2010 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <ostream>
#include <fstream>
#include <utility>

#include <dawg/utils.h>
#include <boost/algorithm/string.hpp>

#include <dawg/residue.h>

namespace dawg {

class output {
public:
	output() : do_op(&output::print_aln), p_out(nullptr), format_id(0),
		rep(0), label_width(0),
		do_append(false), do_split(false), do_label(false),
		split_id_offset(0) { }

	bool open(const char *file_name, unsigned int max_rep=0,
		bool split = false, bool append=false, bool label=false);

	inline bool operator()(const alignment& aln) {
		open_next();
		if(p_out == nullptr)
			return false;
		std::ostream &out = *p_out;
		if(!do_split)
			out << ((rep == 0) ? block_head : block_between);
		out << block_before;
		(this->*do_op)(aln);
		out << block_after;
		if(!do_split == 0 && rep == last_rep)
			out << block_tail;
		++rep;
		out.flush();
		return true;
	}

	template<class T>
	bool set_format(T format);

	inline void set_ostream(std::ostream &os) {
		p_out = &os;
	}
	inline void set_ostream(std::ostream *os) {
		p_out = os;
	}

	inline void set_blocks(const char *h, const char *w,
		const char *t, const char *b, const char *a) {
		block_head.assign(h);
		block_between.assign(w);
		block_tail.assign(t);
		block_before.assign(b);
		block_after.assign(a);
	}

private:
	void (output::*do_op)(const alignment& aln);

	void print_aln(const alignment& aln);
	void print_poo(const alignment& aln);
	void print_fasta(const alignment& aln);
	void print_nexus(const alignment& aln);
	void print_phylip(const alignment& aln);

protected:
	bool open_file(const char* file_name);
	bool open_next();

	std::ostream *p_out;
	std::ofstream fout;
	std::size_t format_id;
	unsigned int rep, label_width, last_rep;
	bool do_append, do_split, do_label;
	std::string current_label, split_file_name;
	std::string::size_type split_id_offset;

	std::string block_head;
	std::string block_between;
	std::string block_tail;
	std::string block_before;
	std::string block_after;
};

template<class T>
bool output::set_format(T format) {
	static constexpr char format_keys[][10] = {
		"aln", "poo", "fasta", "fsa",
		"nexus", "phylip"
	};
	static void (output::*format_ops[])(const alignment& aln) = {
		&output::print_aln, &output::print_poo,
		&output::print_fasta, &output::print_fasta,
		&output::print_nexus, &output::print_phylip
	};
	format_id = key_switch(format, format_keys);
	if(format_id == (std::size_t)-1) {
		format_id = 0;
		do_op = &output::print_aln;
		return false;
	}
	do_op = format_ops[format_id];
	return true;
}

} // namespace dawg

#endif // DAWG_OUTPUT_H
