#pragma once
#ifndef DAWG_RESIDUE_H
#define DAWG_RESIDUE_H
/****************************************************************************
 *  Copyright (C) 2009-2018 Reed A. Cartwright, PhD <reed@scit.us>          *
 ****************************************************************************/

#ifndef __STDC_CONSTANT_MACROS
#	define __STDC_CONSTANT_MACROS 1
#endif
#ifndef __STDC_LIMIT_MACROS
#	define __STDC_LIMIT_MACROS 1
#endif

#include <dawg/log.h>
#include <dawg/error.h>

#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cstring>
#include <unordered_map>

#include <boost/cstdint.hpp>
#include <boost/bind.hpp>
#include <boost/range/sub_range.hpp>

namespace dawg {

namespace details {

///////////////////////////////////////////////////////////
/// \brief Aligned Sequence
///////////////////////////////////////////////////////////
struct aligned_sequence {
	std::string label;
	std::string seq;
};

} // details

class residue;
typedef std::vector<residue> sequence;

struct alignment : public std::vector<details::aligned_sequence> {
	std::string::size_type max_label_width;
	int seq_type;
};

///////////////////////////////////////////////////////////
/// \brief Residue model
///	\brief Create function
/// \brief Sequence type
/// \brief Creating a root
/// \brief Random residues
/// \brief Evolve residues
/// \brief Constant or variable residue widths
/// \brief Indel widths as well?
///////////////////////////////////////////////////////////
class residue {
public:
	typedef boost::uint64_t data_type;

	static constexpr data_type base_mask      =  UINT64_C(0x00000000000000FF);
	static constexpr data_type branch_mask    =  UINT64_C(0x0000FFFFFFFFFF00);
	static constexpr data_type rate_mask	  =  UINT64_C(0xFFFF000000000000);
	static constexpr data_type rate_shift     =  48;
	static constexpr data_type base_bit_width =  8;
	static constexpr data_type branch_inc     =  0x100;

	inline data_type base() const { return data_ & base_mask; }
	inline void base(data_type b) { data_ = (b & base_mask) | (data_ & ~base_mask); }

	inline data_type branch() const { return data_ & branch_mask; }
	inline void branch(data_type u) { data_ = (u & branch_mask) | (data_ & ~branch_mask); }

	inline data_type rate_cat() const { return data_ >> rate_shift; }
	inline void rate_cat(data_type k) { data_ = (k << rate_shift) | (data_ & ~rate_mask); }

	inline data_type data()  const { return data_; }
	inline void data(data_type d) { data_ = d; }
	inline void data(data_type n, data_type r, data_type b) {
		data_ = (n & base_mask) | (b & branch_mask) | (r << rate_shift);
	}

	inline bool is_base(data_type u) const { return (base() == (u & base_mask)); }
	inline bool is_branch(data_type u) const { return (branch() == (u & branch_mask)); }

	residue() : data_(0)  { }
	residue(data_type n, data_type r, data_type b) :
		data_((n & base_mask) | (b & branch_mask) | (r << rate_shift))
	{

	}


protected:
	data_type data_;
};

// template<class CharType, class CharTrait>
// std::basic_ostream<CharType, CharTrait>&
// operator<<(std::basic_ostream<CharType, CharTrait>& o, const dawg::residue &v) {
// 	if(!o.good()) return o;
// 	o << v.base();
// 	return o;
// }

///////////////////////////////////////////////////////////
/// \brief Residue Exchange
/// Defines the model and coding functions for the residue
///////////////////////////////////////////////////////////
class residue_exchange {
public:
	enum { MODDNA = 0, MODRNA=2, MODAA=4, MODCOD=6, MODEND=30};
	enum { DNA = 0, AA = 1, CODON = 2};

	typedef residue_exchange self_type;

	typedef boost::sub_range< const char [64] > str_type;

	explicit residue_exchange(int m=DNA) { model(m,0,0,false,false,false); }

	inline bool model(unsigned int type, unsigned int code, bool rna,
			bool lowercase, bool markins, bool keepempty) {

		static constexpr char sIns[] = "-+";
		// table for going from base->char
		// TODO: Allow codons to be translated into aa
		// table is 64 wide and 30 down
		static constexpr char mods[] =
			"ACGT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // DNA
			"acgt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // dna
			"ACGU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // RNA
			"acgu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // rna
			"ACDEFGHIKLMNPQRSTVWY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // AA
			"acdefghiklmnpqrstvwy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-" // aa
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" // codons
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqr!!uvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ@=KL-MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ@-KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-=KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEFGHIJ-=KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"ABCDEFGHIJ-!KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"ABCDEF-HIJ!=KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
			"AB-DEFGHIJ!!KL!MNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
		;
		// tables for going from char->base
		// 1 genetic code is 80 elements long
		static constexpr char rmods[] = {
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1, // DNA & RNA
			-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3, 3,-1,-1,
			-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,
			-1,-1,-1,-1,-1,-1,-1,-1, 3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1, // AA
			 2, 3, 4, 5, 6, 7,-1, 8, 9,10,11,-1,12,13,14,15,16,-1,17,18,
			-1,19,-1,-1,-1,-1,-1,-1,-1, 0,-1, 1, 2, 3, 4, 5, 6, 7,-1, 8,
			 9,10,11,-1,12,13,14,15,16,-1,17,18,-1,19,-1,-1,-1,-1,-1,-1,
			54,55,56,57,58,59,60,61,62,63,-1,-1,-1,11,-1,-1,10, 0, 1, 2, // CODON
			 3, 4, 5, 6, 7, 8, 9,12,13,15,16,14,17,18,19,20,21,22,23,24,
			25,26,27,-1,-1,-1,-1,-1,-1,28,29,30,31,32,33,34,35,36,37,38,
			39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,-1,-1,-1,-1,-1
		};

		markins_ = markins;
		keepempty_ = keepempty;
		lowercase_ = lowercase;
		rna_ = rna;
		type_ = type; // set sequence type DNA, AA, or CODON
		nuc_ = ((rna) ? MODRNA : MODDNA) | ((lowercase) ? 1 : 0);
		code_ = code;

		if(MODCOD+code_ >= MODEND || mods[(MODCOD+code_)*64] == '!') {
			throw dawg::dawg_error_t(dawg_error::invalid_genetic_code);
		}

		switch(type_) {
			case DNA:
				cs_decode_ = &mods[nuc_*64];
				break;
			case AA:
				cs_decode_ = &mods[(MODAA + (lowercase ? 1 : 0))*64];
				break;
			case CODON:
				cs_decode_ = &mods[(MODCOD+code_)*64];
				break;
			default:
				throw dawg::dawg_error_t(dawg_error::invalid_sequence_type);
		}
		cs_encode_ = &rmods[type_*80];
		cs_ins_ = &sIns[(markins_ ? 1 : 0)];

		gap_ = static_cast<unsigned int>(strchr(cs_decode_, '-')-cs_decode_);

		do_op_append =  (type_ == CODON) ?
			&residue_exchange::do_op_append_cod :
			&residue_exchange::do_op_append_res ;
		do_op_appendi = (type_ == CODON) ?
			&residue_exchange::do_op_appendi_cod :
			&residue_exchange::do_op_appendi_res ;

		return true;
	};
	inline bool is_same_type(unsigned int type, bool markins, bool keepempty) const {
		return (type == type_ && markins == markins_ && keepempty == keepempty_);
	}
	inline bool is_same_model(unsigned int type, unsigned int code, bool rna,
			bool lowercase, bool markins, bool keepempty) const {
		return ( type == type_       && code == code_
			  && rna  == rna_        && lowercase == lowercase_
			  && markins == markins_ && keepempty == keepempty_
			  );
	}
	inline bool is_keep_empty() const { return keepempty_; }

	inline residue::data_type gap_base() const { return gap_; }

	inline residue::data_type encode(char ch) const {
		char ret = (ch >= '0') ? (cs_encode_[ch - '0']) : -1;
		return static_cast<residue::data_type>(ret);
	}

	///////////////////////////////////////////////////////////
	/// \brief encode
	/// \param root_seq The root sequence to encode
	///////////////////////////////////////////////////////////
	inline sequence encode(const std::string &root_seq) const {
		static const auto getCodonNumber = []
			(const char a, const char b, const char c)->unsigned int {
				return a | (b << 8) | (c << 16);
			};
		sequence residues;
		if (type_ != CODON) {
			for (size_t i = 0; i != root_seq.size(); ++i) {
				auto base = encode(root_seq.at(i));
				if (base == static_cast<decltype(base)>(-1)) {
					throw dawg::dawg_error_t(dawg_error::invalid_user_sequence);
					return {};
				}
				residues.emplace_back(base, 0, 0);
			}
		} else {
			if (root_seq.size() % 3 != 0) {
				throw dawg::dawg_error_t(dawg_error::invalid_user_sequence, \
				    std::string("sequence does not fit codon"));
				return {};
			}
			for (size_t i = 0; i + 2 < root_seq.size(); i += 3) {
				residues.emplace_back(triplet_to_codon(
					getCodonNumber(root_seq.at(i), root_seq.at(i + 1), root_seq.at(i + 2))), 0, 0);
				if (residues.back().data() == -1) {
					throw dawg::dawg_error_t(dawg_error::invalid_user_sequence);
					return {};
				}
			}
		}
		return residues;
	}

	inline char decode(residue::data_type r) const {
		return cs_decode_[r & 63];
	}

	inline char decode(const residue &r) const {
		return decode(r.base());
	}

	inline char decode_ins() const {
		return cs_ins_[0];
	}

	/// \brief codon number -> cod64
	static inline char codon_to_cod64(unsigned int p) {
		const char s[] = "ABCDEFGHIJ@=KLOMNPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
		return s[p&63];
	}

	/// \brief cod64 -> codon number
	static inline unsigned int cod64_to_codon(char c) {
		static constexpr char a[] = {
			// cod64 -> codon number
			54,55,56,57,58,59,60,61,62,63,-1,-1,-1,11,-1,-1,10, 0, 1, 2,
			 3, 4, 5, 6, 7, 8, 9,12,13,15,16,14,17,18,19,20,21,22,23,24,
			25,26,27,-1,-1,-1,-1,-1,-1,28,29,30,31,32,33,34,35,36,37,38,
			39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,-1,-1,-1,-1,-1
		};
		return (c >= '0') ? a[c-'0'] : -1;
	}

	/// \brief codon number -> triplet
	static inline unsigned int codon_to_triplet(unsigned int p, int type=0) {
		const char ss[] = "TCAGtcagUCAGucag";
		const char *s = &ss[(type&3)*4];
		unsigned int u = s[p&3];
		u = (u << 8) | s[(p>>2)&3];
		u = (u << 8) | s[(p>>4)&3];
		return u;
	}

	/// \brief triplet -> codon number
	static inline unsigned int triplet_to_codon(unsigned int p) {
		// randomly cryptic code:
		// In ncbi descriptions of genetic code,
		//     T|U=0,C=1,A=2,G=3
		// Let f(p) = (p&6)/2.  For ASCII:
		//     f(T|U)= 2, f(C) = 1, f(A)=0, f(G)=3
		// To convert ascii to ncbi: g(p)= (3*p+2)%4
		unsigned int u = p & 0x060606;
		u = (u+(u/2)+0x020202) & 0x030303;
		// reverse bytes and pack into lower 8 bits
		// replaces three shifts, three masks, and two ors
		return ((u*1049601)/65536) % 64;
	}

	void append_residue(std::string &ss, const residue &r) const {
		return (this->*do_op_append)(ss,r);
	}

	void append_ins(std::string &ss) const {
		return (this->*do_op_appendi)(ss);
	}

	/// \brief get the protein code from a genetic code
	inline static const char* get_protein_code(unsigned int code) {
		static constexpr char s[] =
			"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
			"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
			"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
		;
		return &s[64*(code%24)];
	}

protected:
	void (residue_exchange::*do_op_append)(std::string &ss, const residue &r) const;
	void (residue_exchange::*do_op_appendi)(std::string &ss) const;

	void do_op_append_res(std::string &ss, const residue &r) const {
		ss.append(1, decode(r));
	}
	void do_op_append_cod(std::string &ss, const residue &r) const {
		unsigned int b = r.base()&63;
		if(b == gap_)
			ss.append(3, '-');
		else if(cs_decode_[b] == '!')
			ss.append(3, '*');
		else {
			unsigned int u = codon_to_triplet(b, nuc_);
			ss.append(1, (char)(u));
			ss.append(1, (char)(u>>8));
			ss.append(1, (char)(u>>16));
		}
	}
	void do_op_appendi_res(std::string &ss) const {
		ss.append(1, decode_ins());
	}
	void do_op_appendi_cod(std::string &ss) const {
		ss.append(3, decode_ins());
	}

	unsigned int type_, code_, nuc_, gap_;
	bool rna_, markins_, keepempty_, lowercase_;
	const char *cs_decode_, *cs_ins_, *cs_encode_;
};

} // namespace dawg
#endif
