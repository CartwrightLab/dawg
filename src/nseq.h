#pragma once
#ifndef DAWG_NSEQ_H
#define DAWG_NSEQ_H

#include <vector>
#include <list>
#include <algorithm>
#include <functional>
#include <iterator>

namespace dawg {
// http://www.cs.princeton.edu/courses/archive/fall08/cos226/lectures/10BalancedTrees-2x2.pdf

template<class _T, class _W> class finger_tree;

namespace detail {

template<class _T, class _W>
class finger_tree_node_iterator
	: public std::iterator< std::bidirectional_iterator_tag,
	                        typename dawg::finger_tree<_T, _W>::node >
{
public:
	typedef finger_tree_node_iterator<_T, _W> self_type;
	typedef typename dawg::finger_tree<_T, _W>::node node_type;
	typedef std::iterator<std::bidirectional_iterator_tag, node_type> base_type;
	
	typedef typename base_type::value_type value_type;
	typedef typename base_type::difference_type difference_type;
	typedef typename base_type::pointer pointer;
	typedef typename base_type::reference reference;
	typedef typename base_type::iterator_category iterator_category;
	
	finger_tree_node_iterator() : p_node(node_type::null) { }
	explicit finger_tree_node_iterator(pointer p) : p_node(p) { }
	
	reference operator*() const { return *p_node; }
	pointer operator->() const { return p_node;	}
	self_type& operator++() {
		if(p_node->right != NULL) {
			p_node = p_node->right;
			while(p_node->left != 0)
				p_node = p_node->left;
		} else {
			pointer p_up = p_node->up;
			while( p_up != NULL && p_node == p_up->right) {
				p_node = p_up;
				p_up = p_node->up;
			}
			
		}
		return *p_node;
	}
	self_type& operator++(int) {
		
	}
	self_type& operator--() {
		
	}
	self_type& operator--(int) {
		
	}
	bool operator==(const self_type& x) const {
		return p_node == x.p_node;
	}
	bool operator!=(const self_type& x) const {
		return p_node != x.p_node;
	}
	
protected:
	pointer p_node;
};

}; // namespace dawg::detail

template<class _T, class _W>
class finger_tree {
public:
	typedef finger_tree<_T, _W> self_type;
	struct node;
	typedef _T data_type;
	typedef _W weight_type;
	typedef detail::finger_tree_node_iterator<_T, _W> iterator;

	finger_tree() {
		head.up = &head;
		head.left = &head;
		head.right = &head;
	}
	~finger_tree() {
		if(head.up != &head)
			delete head.up;
		head.left = NULL;
		head.right = NULL;
	}
	
	iterator begin() {
		return iterator(head.left);
	}
	iterator end() {
		return iterator(&head);
	}
	iterator root() {
		return iterator(head.up);
	}
	
	struct node {
		typedef node* pointer;
		
		pointer left, right, up;
		bool color; // red == true; black == false
		data_type val;
		weight_type weight;

		node() : left(NULL), right(NULL), up(NULL), color(false) { }
		node(const data_type &v) : left(NULL), right(NULL), up(NULL), val(v), color(true) { }
		~node() {
			if(left != NULL)
				delete left;
			if(right != NULL)
				delete right;
		}

		pointer rotate_left() {
			pointer x = right;
			right = x->left;
			x->left = this;
			x->color = color;
			color = true;
			if(up->up == this)
				up->up = x;
			else if(up->left == this)
				up->left = x;
			else
				up->right = x;
			x->up = up;
			up = x;
			return x;
		}
	
		pointer rotate_right() {
			pointer x = left;
			left = x->right;
			x->right = this;
			x->color = color;
			color = true;
			if(up->up == this)
				up->up = x;
			else if(up->left == this)
				up->left = x;
			else
				up->right = x;
			x->up = up;
			up = x;
			return x;
		}
	
		void flip_colors() {
			color = true;
			left->color = false;
			right->color = false;
		}

	};
	
	// insert value before position in tree
	iterator insert(iterator pos, const data_type &val) {
		typename node::pointer p = &(*pos);
		typename node::pointer x = new node(val);

		// find null link that is just to the left of p
		// update start and end links if needed
		// check of empty tree
		if(begin() == end()) {
			head.left = x;
			head.right = x;
			head.up = x;
			x->color = false;
		} else if(is_null(p->left)) {
			p->left = x;
			if(p == head.left)
				head.left = x;
		} else {
			if(pos == end())
				p = head.right;
			else for(p = p->left; !is_null(p->right); p = p->right)
				/*noop*/;
			p->right = x;
			if(p == head.right)
				head.right = x;
		}
		x->up = p;

		for(; p != &head; p = p->up) {
			if(is_red(p->right) && !is_red(p->left))
				p = p->rotate_left();
			if(is_red(p->left) && is_red(p->left->left))
				p = p->rotate_right();
			if(is_red(p->left) && is_red(p->right))
				p->flip_colors();
		}
		head.up->color = false;
		
		return iterator(x);
	}
			
	
private:
	inline bool is_null(typename node::pointer p) const {
		return (p == NULL);
	}
	inline bool is_red(typename node::pointer p) const {
		return (!is_null(p) && p->color == true);
	}
	
	node head;
};

};

class ResidueFactory
{
public:
	typedef unsigned int model_type;
	static const model_type modelDNA = 0;
	static const model_type modelRNA = 1;
	static const model_type modelAA  = 2;
	static const model_type modelCod = 3;

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
	typedef Block::const_iterator const_seq_iterator;
	
	template<typename It>
	seq_iterator operator()(It b, It e) {
		// should speed test
		store.push_back(Block());
		store.back().reserve(e-b);
		std::transform(b, e, std::back_inserter(store.back()),
			std::bind1st(std::mem_fun(&ResidueFactory::makeResidue), this));
		return store.back().begin();
	}
	seq_iterator operator()(seq_iterator b, seq_iterator e) {
		store.push_back(Block(b,e));
		return store.back().begin();
	}
	seq_iterator operator()(const_seq_iterator b, const_seq_iterator e) {
		store.push_back(Block(b,e));
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
		return Residue(encode(ch), 1.0);
	}
	Residue makeResidue(char ch, double d) const {
		return Residue(encode(ch), d);
	}
	
	Residue::base_type encode(char ch) const {
		static Residue::base_type dna[] = {0,1,3,2};
		if(seq_model == modelDNA || seq_model == modelRNA)
			return dna[(ch&6u) >> 1];
		return ~0;
	}
	char decode(Residue::base_type r) const {
		static char dna[] = "ACGT";
		static char rna[] = "ACGU";
		// Include O & U, two rare amino acids
		static char aa[] = "ACDEFGHIKLMNPQRSTVWYOU";
		// Encode 64 possible codons using a base-64 notation
		// b/c this function returns a single char.
		// These can later be converted to three nucleotides or 1 amino acid
		static char cod[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.";

		switch(seq_model) {
		case modelDNA:
			return dna[r];
		case modelRNA:
			return rna[r];
		case modelAA:
			return aa[r];
		case modelCod:
			return cod[r];
		default:
			break;
		}
		return '?';
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
	residue_iterator(high_iterator h, high_iterator e)
		: hit(h), hit_end(e) {
		while(dirtyhit())
			++hit;
		lit = (hit != hit_end) ? hit->begin() : low_iterator();
	}
	residue_iterator(high_iterator h, high_iterator(e), low_iterator l)
		: hit(h), hit_end(e), lit(l) {
		/* assumes consistency */
	}

	inline bool dirtyhit() {
		return ( hit != hit_end && (hit->empty() || hit->IsDeletion() ));
	}

	this_type& operator++() {
		if(++lit == hit->end()) {
			do { ++hit; } while(dirtyhit());
			lit = (hit != hit_end) ? hit->begin() : low_iterator();
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
			do {--hit; } while(dirtyhit());
			lit = (hit != hit_end) ? --hit->end() : low_iterator();
		}
		return *this;
	}
	this_type operator--(int) {
		this_type ret = *this;
		--(*this);
		return ret;
	}
	bool operator==(const this_type &other) const {
		return (hit==other.hit && (hit == hit_end || lit == other.lit));
	}
	bool operator!=(const this_type &other) const {
		return !(*this == other);
	}

	reference operator*() const {
		return *lit;
	}
	pointer operator->() const {
		return (&**this);
	}

	high_iterator _hit() { return hit; }
	low_iterator _lit() { return lit; }
	
protected:
	high_iterator hit;
	high_iterator hit_end;
	low_iterator lit;
};

class Sequence2 {
public:
	// A sequence contains two important pieces of data
	// 1. a vector of vectors of nucleotides.  This contains the raw Part data
	// 2. a vector of Part information
	
	static ResidueFactory makeSeq;
		
	class Part {
	public:
		typedef unsigned int color_type;
		typedef ResidueFactory::Block data_type;
		typedef data_type::value_type value_type;
		typedef data_type::iterator iterator;
		typedef data_type::const_iterator const_iterator;
		typedef data_type::reference reference;
		typedef data_type::const_reference const_reference;
		typedef data_type::size_type size_type;
		
		static const color_type BranchInc   =  0x2;
		static const color_type MaskBranch	= ~0x1;
		static const color_type MaskType	=  0x1; // 1
		static const color_type TypeDel		=  0x1; // 1
		static const color_type TypeExt		=  0x0; // 0
		
		Part(iterator it, size_type len, color_type uT, color_type uB) :
			itBegin(it), szLength(len), color((uT & MaskType) | (uB & MaskBranch))
			{ }
		Part(iterator it, size_type len, color_type u) :
			itBegin(it), szLength(len), color(u)
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
		
		inline bool empty() const { return (szLength == 0); }

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
		
		inline Part split(iterator loc) {
			Part ret;
			if(loc == itBegin)
				return ret;
			ret.itBegin = loc;
			ret.szLength = itBegin+szLength-loc;
			ret.color = color;
			szLength -= ret.szLength;
			return ret;			
		}
		
		inline Part clone() const {
			return Part(makeSeq(begin(), end()), szLength, color);
		}
			
	protected:
		iterator itBegin;
		size_type szLength;
		color_type color;

	private:
		Part() : szLength(0) { };
		//Part(const Part &p);
		//Part & operator = (const Part &p);
		
		friend class Sequence2;
	};

	typedef std::list<Part> Parts;
	typedef residue_iterator<Parts::iterator, Part::iterator> iterator;
	typedef residue_iterator<Parts::const_iterator, Part::const_iterator> const_iterator;
	typedef Part::size_type size_type;

	Sequence2() : szLength(0) {
	}

	template<typename It1, typename It2>
	void append(It1 b, It1 e, It2 r) {
		ResidueFactory::seq_iterator it = makeSeq(b,e,r);
		size_type len = e-b;
		parts.push_back(Part(it, len, Part::TypeExt, 0));
		szLength += len;
	}

	template<typename It1>
	void append(It1 b, It1 e) {
		ResidueFactory::seq_iterator it = makeSeq(b,e);
		size_type len = e-b;
		parts.push_back(Part(it, len, Part::TypeExt, 0));
		szLength += len;
	}
	
	void append(const Sequence2 &s2) {
		parts.insert(parts.end(), s2.parts.begin(), s2.parts.end());
		szLength += s2.size();
	}
	
	void insert(size_type pos, const Sequence2 &s2) {
		iterator it = mid(pos);
		Parts::iterator hit = it._hit();
		if(hit == parts.end()) {
			parts.insert(hit, s2.parts.begin(), s2.parts.end());
		} else {
			Part tail = hit->split(it._lit());
			parts.insert(++hit, s2.parts.begin(), s2.parts.end());
			if(!tail.empty())
				parts.insert(hit, tail);
		}
		szLength += s2.size();
	}
	
	void erase(size_type pos, size_type len) {
		iterator it = mid(pos);
		Parts::iterator hit = it._hit();
		if(hit == parts.end())
			return; //nothing to do
		// find start location and split if needed
		Part tail = hit->split(it._lit());
		if(!tail.empty()) {
			parts.insert(++hit, tail);
			--hit;
		}
		// hit now contains the first part of the deletion
		// find the last part
		while(hit != parts.end() && hit->size() <= len) {
			// mark all of this part as deleted
			hit->SetType(Part::TypeDel);
			len -= hit->size();
			// advance past deletions
			do { ++hit; } while(hit->IsDeletion());
		}
		if(hit == parts.end() || len == 0)
			return; // done
		// split again
		Part tail2 = hit->split(hit->begin()+len);
		hit->SetType(Part::TypeDel);
		parts.insert(++hit, tail2);
	}
	
	void optimize() {
		if(parts.empty())
			return;
		Parts::iterator itpp = parts.begin();
		Parts::iterator it = itpp++;
		while(itpp != parts.end()) {
			if(it->end() == itpp->begin() && it->GetColor() == itpp->GetColor()) {
				it->szLength += itpp->szLength;
				itpp = parts.erase(itpp);
			} else {
				++it;
				++itpp;
			}
		}
	}
	
	Sequence2 clone() {
		Sequence2 ret;
		ret.szLength = szLength;
		std::transform(parts.begin(), parts.end(), std::back_inserter(ret.parts),
			std::mem_fun_ref(&Part::clone));
		return ret;
	}

	size_type size() const { return szLength; }

	// begin, end, and mid cycle through extant residues

	inline iterator begin() {
		return iterator(parts.begin(),parts.end());
	}
	inline iterator end() {
		return iterator(parts.end(),parts.end());
	}
	inline const_iterator begin() const {
		return const_iterator(parts.begin(), parts.end());
	}
	inline const_iterator end() const {
		return const_iterator(parts.end(),  parts.end());
	}
		
	inline iterator mid(size_type pos) {
		if(pos >= size())
			return end();
		Parts::iterator it = parts.begin();
		for(; it != parts.end(); ++it) {
			if(it->IsDeletion())
				continue;
			if(pos < it->size())
				break;
			 pos -= it->size();
		}
		return (it == parts.end()) ? iterator(parts.end(), parts.end())
			: iterator(it, parts.end(), it->begin()+pos);
	}
	
	const Parts& _parts() const { return parts; }

protected:
	Parts parts;
	size_type szLength;
};


#endif
