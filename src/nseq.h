#pragma once
#ifndef DAWG_NSEQ_H
#define DAWG_NSEQ_H

#include <vector>
#include <list>
#include <algorithm>
#include <functional>

class ResidueFactory
{
public:
	typedef unsigned int model_type;
	static const model_type modelDNA = 0;
	static const model_type modelRNA = 1;
	
	struct Residue {
		typedef unsigned char base_type;
		typedef float rate_type;
		
		base_type base;
		rate_type rate;
		
		Residue() { }
		
		Residue(base_type b, rate_type r)
			: base(b), rate(r) { }
	};

	typedef std::vector<Residue> Block; 
	typedef std::list<Block> DataStore;
	typedef Block::iterator seq_iterator;
	
	template<typename It>
	seq_iterator operator()(It b, It e) {
		// should speed test
		store.push_back(Block());
		store.back().reserve(e-b);
		std::transform(b, e, std::back_inserter(store.back()),
			std::bind1st(std::mem_fun(&ResidueFactory::makeResidue), this));
		return store.back().begin();
	}

	struct biop : public std::binary_function<char, double, Residue> {
		biop(const ResidueFactory *o) : obj(o) { }
		const ResidueFactory *obj;
		Residue operator()(char ch, double d) {
			return obj->makeResidue(ch,d);
		}
	};
	template<typename It1, typename It2>
	seq_iterator operator()(It1 b, It1 e, It2 r) {
		store.push_back(Block());
		store.back().reserve(e-b);
		std::transform(b, e, r, std::back_inserter(store.back()), biop(this));
		return store.back().begin();
	}

	inline void SetModel(model_type a) { seq_model = a; };
	
	Residue makeResidue(char ch) const {
		return Residue(translate(ch), 1.0);
	}
	Residue makeResidue(char ch, double d) const {
		return Residue(translate(ch), d);
	}	
	
	Residue::base_type translate(char ch) const {
		static Residue::base_type dna[] = {0,1,3,2};
		
		if(seq_model == modelDNA || seq_model == modelRNA)
			return dna[(ch&6u) >> 1];
		return ~0;
	}
	
	
	void clear() { store.clear(); }

protected:
	DataStore store;
	model_type seq_model;
};

template<typename ith, typename itl>
class residue_iterator {
public:
	typedef ith high_iterator;
	typedef itl low_iterator;
	typedef typename low_iterator::reference reference;
	typedef typename low_iterator::pointer pointer;
	typedef typename low_iterator::value_type value_type;
	typedef std::bidirectional_iterator_tag iterator_category;

	typedef residue_iterator<high_iterator, low_iterator> this_type;
	
	residue_iterator() {}
	residue_iterator(high_iterator h, bool b=true)
		: hit(h), ignDel(b) {
		for(;ignDel && hit->IsDeletion();++hit)
			/*noop*/;
		lit = hit->begin();
	}
	residue_iterator(high_iterator h, low_iterator l, bool b=true)
		: hit(h), lit(l), ignDel(b) {
		/* assumes consistency */
	}

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
			do {--hit; } while(ignDel && hit->IsDeletion());
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
	bool ignDel;
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
		
		static const color_type BranchInc   =  0x2;
		static const color_type MaskBranch	= ~0x1;
		static const color_type MaskType	=  0x1; // 1
		static const color_type TypeDel		=  0x1; // 1
		static const color_type TypeExt		=  0x0; // 0
		
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
		inline bool IsDeletion() const { return IsType(TypeDel); }
		inline bool IsExtant()   const { return IsType(TypeExt); }
		
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
	typedef residue_iterator<Parts::iterator, Part::iterator> iterator;
	typedef residue_iterator<Parts::const_iterator, Part::const_iterator> const_iterator;
	typedef Part::size_type size_type;

	Sequence2() : szLength(0), szAlnLength(0) {
		// Push back termination sequence
		parts.push_back(Part(Part::iterator(), 0, 0, 0));
	}

	size_type size() const { return szLength; }
	size_type asize() const { return szAlnLength; }

	// begin, end, and mid cycle through extant residues

	inline iterator begin(bool b=true) {
		return iterator(parts.begin(),b);
	}
	inline iterator end(bool b=true) {
		return iterator(--parts.end(), b);
	}
	inline iterator mid(size_type pos, bool b=true) {
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
		return iterator(it, it->begin()+pos,b);
	}

	inline const_iterator begin(bool b=true) const {
		return const_iterator(parts.begin(),b);
	}
	inline const_iterator end(bool b=true) const {
		return const_iterator(--parts.end(),b);
	}
	inline const_iterator mid(size_type pos, bool b=true) const {
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
		return const_iterator(it, it->begin()+pos,b);
	}

protected:
	Parts parts;
	size_type szLength;
	size_type szAlnLength;
};

#endif
