#pragma once
#ifndef DAWG_NSEQ_H
#define DAWG_NSEQ_H

#include <vector>
#include <list>
#include <algorithm>


class ResidueFactory
{
public:
	typedef unsigned int model_type;
	static const model_type modelDNA = 0;

	struct Residue {
		typedef unsigned char base_type;
		typedef float rate_type;
		
		base_type base;
		rate_type rate;
	};

	typedef std::vector<Residue> Block; 
	typedef std::list<Block> DataStore;
	typedef Block::iterator seq_iterator;
	
	template<typename It>
	seq_iterator operator()(It b, It e) {
		// should speed test
		store.push_back(Block(e-b));
		//store.back().resize(e-b);
		stl::transform(b, e, store.back().begin(), );
	}

	inline void SetModel(model_type a) { seq_model = a; };
	
	Residue::base_type translate(char ch) {
		static Residue::base_type dna[] = {0,1,3,2};
		
		if(seq_model == modelDNA)
			return dna[(ch&6u) >> 1]; // ACGT/acgt -> 0132
		return ~0;
	}

protected:
	DataStore store;
	model_type seq_model;
};

template<typename ith, typename itl>
class nested_iterator {
public:
	typedef ith high_iterator;
	typedef itl low_iterator;
	typedef typename low_iterator::reference reference;
	typedef typename low_iterator::pointer pointer;

	typedef nested_iterator<high_iterator, low_iterator> this_type;
	
	nested_iterator() {}
	nested_iterator(high_iterator h, low_iterator l) : hit(h), lit(l) { }

	this_type& operator++() {
		if(++lit == hit->end())
			lit = (++hit)->begin();
		return *this;
	}
	this_type operator++(int) {
		this_type ret = *this;
		++(*this);
		return ret;
	}
	this_type& operator--() {
		if(lit == hit->begin())
			--(lit = (--hit)->end());
		return *this;
	}
	this_type operator--(int) {
		this_type ret = *this;
		--(*this);
		return ret;
	}
	bool operator==(this_type &other) {
		return (hit==other.hit && lit == other.lit);
	}

	reference operator*() const {
		return *lit;
	}
	pointer operator->() const {
		return (&**this);
	}
	
protected:
	high_iterator hit;
	low_iterator lit;
};

template<typename ith, typename itl>
class nested_seq_iterator {
public:
	typedef ith high_iterator;
	typedef itl low_iterator;
	typedef typename low_iterator::reference reference;
	typedef typename low_iterator::pointer pointer;

	typedef nested_seq_iterator<high_iterator, low_iterator> this_type;
	
	nested_seq_iterator() {}
	nested_seq_iterator(high_iterator h, low_iterator l) : hit(h), lit(l) { }

	this_type& operator++() {
		if(++lit == hit->end()) {
			do { ++hit; } while(hit->IsDeletion());
			lit = hit->begin();
		}
		return *this;
	}
	this_type operator++(int) {
		this_type ret = *this;
		++(*this);
		return ret;
	}
	this_type& operator--() {
		if(lit == hit->begin()) {
			do {--hit; } while(hit->IsDeletion());
			--(lit = hit->end());
		}
		return *this;
	}
	this_type operator--(int) {
		this_type ret = *this;
		--(*this);
		return ret;
	}
	bool operator==(this_type &other) {
		return (hit==other.hit && lit == other.lit);
	}

	reference operator*() const {
		return *lit;
	}
	pointer operator->() const {
		return (&**this);
	}
	
protected:
	high_iterator hit;
	low_iterator lit;
};

class Sequence2
{
public:
	// A sequence contains two important pieces of data
	// 1. a vector of vectors of nucleotides.  This contains the raw Part data
	// 2. a vector of Part information
	
	static ResidueFactory makeSeq;
		
	class Part
	{
	public:
		typedef unsigned int color_type;
		typedef ResidueFactory::Block data_type;
		typedef data_type::value_type value_type;
		typedef data_type::iterator iterator;
		typedef data_type::const_iterator const_iterator;
		typedef data_type::size_type size_type;
		
		static const color_type BranchAdd   = 0x4;
		static const color_type MaskBranch	= ~0x3;
		static const color_type MaskType	= 0x3; // 11
		static const color_type MaskDel		= 0x2; // 10
		static const color_type MaskIns		= 0x1; // 01
		static const color_type TypeRoot	= 0x0; // 00
		static const color_type TypeIns		= 0x1; // 01
		static const color_type TypeDel		= 0x2; // 10
		static const color_type TypeDelIns	= 0x3; // 11
		
		Part(iterator it, size_type len, color_type uT, color_type uB) :
			itBegin(it), szLength(len), color((uT & MaskType) | (uB & MaskBranch))
			{ }
		
		inline color_type GetBranch() const { return color & MaskBranch; }
		inline color_type GetType()   const { return color & MaskType; }
		inline color_type GetColor()  const { return color; }
		
		inline void SetType(color_type u)   { color =  (u & MaskType)   | (color & ~MaskType); }
		inline void SetBranch(color_type u) { color =  (u & MaskBranch) | (color & ~MaskBranch); }
		
		inline void SetColor(color_type uT, color_type uB)
			{ color =  (uT & MaskType) | (uB & MaskBranch); }
		inline void SetColor(color_type u) { color = u; }
		inline bool IsType(color_type u) const { return (GetType() == (u & MaskType)); }
		inline bool IsBranch(color_type u) const { return (GetBranch() == (u & MaskBranch)); }
		inline bool IsDeletion() const { return ((color & MaskDel) == MaskDel); }
		inline bool IsInsertion() const { return ((color & MaskIns) == MaskIns); }
		
		inline iterator begin() { return itBegin;}
		inline iterator end() { return itBegin+szLength; }
		inline iterator mid(size_type pos) {
			//assert pos < length ?
			return itBegin+pos;
		}

		inline const_iterator begin() const { return itBegin; }
		inline const_iterator end() const { return itBegin+szLength; }
		inline const_iterator mid(size_type pos) const {
			//assert pos < length ?
			return itBegin+pos;
		}
		
		inline size_type size() const { return szLength; }
		inline iterator shrink(size_type sz) {
			if(sz <= szLength)
				szLength = sz;
			//else throw error?
			return itBegin+sz;
		}
		
	protected:
		iterator itBegin;
		size_type szLength;
		color_type color;

	private:
		//Part(const Part &p);
		Part & operator = (const Part &p);
	};

	typedef std::list<Part> Parts;
	typedef nested_seq_iterator<Parts::iterator, Part::iterator> iterator;
	typedef nested_seq_iterator<Parts::const_iterator, Part::const_iterator> const_iterator;
	typedef nested_iterator<Parts::iterator, Part::iterator> aln_iterator;
	typedef nested_iterator<Parts::const_iterator, Part::const_iterator> const_aln_iterator;
	typedef Part::size_type size_type;

	Sequence2() : szLength(0), szAlnLength(0) {
		// Push back termination sequence
		parts.push_back(Part(Part::iterator(), 0, 0, 0));
	}

	size_type size() const { return szLength; }
	size_type asize() const { return szAlnLength; }

	// begin, end, and mid cycle through extant residues

	inline iterator begin() {
		Parts::iterator first = parts.begin();
		while(first->IsDeletion()) ++first;
		return iterator(first, first->begin());
	}
	inline iterator end() {
		Parts::iterator last = --parts.end();
		return iterator(last, last->begin());
	}
	inline iterator mid(size_type pos) {
		if(pos >= size())
			return end();
		Parts::iterator it = parts.begin();
		for(; /*it != parts.end()*/; ++it) {
			if(it->IsDeletion())
				continue;
			if(pos < it->size())
				break;
			 pos -= it->size();
		}	
		return iterator(it, it->begin()+pos);
	}

	inline const_iterator begin() const {
		return const_iterator(parts.begin(), parts.begin()->begin());
	}
	inline const_iterator end() const {
		Parts::const_iterator last = --parts.end();
		return const_iterator(last, last->begin());
	}
	inline const_iterator mid(size_type pos) const {
		if(pos >= size())
			return end();
		Parts::const_iterator it = parts.begin();
		for(; /*it != parts.end()*/; ++it) {
			if(it->IsDeletion())
				continue;
			if(pos < it->size())
				break;
			 pos -= it->size();
		}	
		return const_iterator(it, it->begin()+pos);
	}

	// abegin, amid, and aend cycle through the alignment
	inline aln_iterator abegin() {
		Parts::iterator first = parts.begin();
		return aln_iterator(first, first->begin());
	}
	inline aln_iterator aend() {
		Parts::iterator last = --parts.end();
		return aln_iterator(last, last->begin());
	}
	inline aln_iterator amid(size_type pos) {
		if(pos >= asize())
			return aend();
		Parts::iterator it = parts.begin();
		for(; pos < it->size() /*&& it != parts.end()*/; ++it) {
			 pos -= it->size();
		}	
		return aln_iterator(it, it->begin()+pos);
	}
	inline const_aln_iterator abegin() const {
		Parts::const_iterator first = parts.begin();
		return const_aln_iterator(first, first->begin());
	}
	inline const_aln_iterator aend() const {
		Parts::const_iterator last = --parts.end();
		return const_aln_iterator(last, last->begin());
	}
	inline const_aln_iterator amid(size_type pos) const {
		if(pos >= asize())
			return aend();
		Parts::const_iterator it = parts.begin();
		for(; pos < it->size() /*&& it != parts.end()*/; ++it) {
			 pos -= it->size();
		}	
		return const_aln_iterator(it, it->begin()+pos);
	}

protected:
	Parts parts;
	size_type szLength;
	size_type szAlnLength;
};

#endif
