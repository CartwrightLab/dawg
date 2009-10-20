// tree.h - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_TREE_H
#define DAWG_TREE_H

#include "dawg.h"
#include "indel.h"
#include "matrix.h"
#include "sequence.h"
#include "bitree.h"
#include "rand.h"

using namespace dawg;

#include <list>
#include <stack>

// A class used to represent a node in a Newick tree
class NewickNode {
public:
	NewickNode(NewickNode* p, const char *cs, double d);
	std::string m_ssLabel;
	double m_dLen;
	std::auto_ptr<NewickNode> m_pSib;
	std::auto_ptr<NewickNode> m_pSub;

protected:
	void MakeName();
};

// A class to represent nucleotides
//class Nucleotide
//{
//public:
//	typedef unsigned short data_type;
//protected:
//	// First two bits specify base
//	// Second two bits specify type
//	data_type   m_ucNuc;
//	float m_dRate; // 0.0 means invarant

//public:
//	Nucleotide() : m_ucNuc(0xF), m_dRate(1.0) { }
//	Nucleotide(data_type nuc, double rate) : m_ucNuc(nuc), m_dRate((float)rate) { }

//	static const data_type MaskColor	= ~0x7;
//	static const data_type MaskBase		=  0x3; // 011
//	static const data_type MaskType		=  0x4; // 100
//	static const data_type TypeDel		=  0x4; // 100
//	static const data_type TypeExt      =  0x0; // 000
//	static const data_type ColorInc     =  0x8;
//
//	inline data_type GetBase()  const  { return m_ucNuc & MaskBase; }
//	inline data_type GetType()  const  { return m_ucNuc & MaskType; }
//	inline data_type GetColor() const  { return m_ucNuc & MaskColor; }
//	inline void SetBase(data_type uc)  { m_ucNuc =  (uc & MaskBase) | (m_ucNuc & ~MaskBase); }
//	inline void SetType(data_type uc)  { m_ucNuc =  (uc & MaskType) | (m_ucNuc & ~MaskType); }
//	inline void SetColor(data_type uc) { m_ucNuc =  (uc & MaskColor) | (m_ucNuc & ~MaskColor); }
//	inline void SetNuc(data_type ucB, data_type ucT, data_type ucC)
//		{ m_ucNuc =  (ucB & MaskBase) | (ucT & MaskType) | (ucC & MaskColor); }
//	inline void SetNuc(data_type uc) { m_ucNuc = uc; }
//	inline bool IsType(data_type uc) const { return (GetType() == uc); }
//	inline bool IsDeleted() const { return (GetType() == TypeDel); }
//	inline bool IsExtant()  const { return (GetType() == TypeExt); }

//	inline double GetRate() const { return m_dRate; }
//	inline void SetRate(double r) { m_dRate = (float)r; }
//
//	bool FromChar(char ch);
//	char ToChar() const;

//};

// A class that represent a sequence of nucleotides
//class Sequence : public std::vector<Nucleotide>
//{
//public:
//	typedef std::vector<Nucleotide> Base;
//	Sequence() : m_uLength(0) { }
//	explicit Sequence(size_type uSize) : Base(uSize, Nucleotide(0xF, -1.0))
//	{
//		m_uLength = uSize;
//	}
//	size_type SeqLength() const { return m_uLength; }
//
//	// find the uPos-th true nucleotide (skips gaps)
//	const_iterator SeqPos(size_type uPos) const;
//	iterator SeqPos(size_type uPos);

//	size_type Insertion(iterator itPos, const_iterator itBegin, const_iterator itEnd);
//	size_type Deletion(iterator itBegin, Base::size_type uSize);

//	void Append(const Sequence &seq);

//	void ToString(std::string &ss) const;

//private:
//	size_type m_uLength;
//};

// The recombinant tree data structure
struct AlignData;

class Tree
{
public:

	// A node in the tree
	struct Node
	{
		typedef std::vector<dawg::residue> Sequence;
		typedef Sequence::value_type Nucleotide;
		typedef std::map<std::string, Tree::Node> Map;
		typedef Sequence::size_type size_type;
		Sequence m_vSeq;
		std::vector<std::string> m_vAncestors;
		double m_dBranchLen;
		std::string m_ssName;
		bool m_bTouched;

		Node() : m_bTouched(false), m_dBranchLen(0.0) { }
	};
	typedef Node::Sequence Sequence;
	typedef Node::Nucleotide Nucleotide;
	typedef Sequence::size_type size_type;

	typedef std::map<std::string, std::string> Alignment;

	// Setup the model of evolution
	bool SetupEvolution(double pFreqs[], double pSubs[],
		const IndelModel::Params& rIns, const IndelModel::Params& rDel,
		double dGamma, double dIota, double dTreeScale,
		int uKeepFlank);

	// Setup the root node
	bool SetupRoot(const std::vector<std::string> &vSeqs, const std::vector<unsigned int> &vData,
		const std::vector<std::vector<double> > &vRates);

	// Draw a random relative rate of substitution from the evolutionary parameters
	inline Nucleotide::rate_type RandomRate() const {
		if(AreRatesConstant())
			return Nucleotide::rate_type(1.0);
		if(m_dIota > DBL_EPSILON && rand_bool(m_dIota))
			return Nucleotide::rate_type(0.0);  // Site Invariant
		return static_cast<Nucleotide::rate_type>(rand_gamma1(m_dGamma)); // Gamma with mean 1.0 and var of m_dGamma
	}
	bool AreRatesConstant() const {
		return (m_dIota < DBL_EPSILON && m_dGamma < DBL_EPSILON);
	}
	// Draw a random base from the evolutionary parameters
	inline Nucleotide::data_type RandomBase() const {
		return m_dFreqsCum(rand_real());
	}

	// Draw a random nucleotide (base and rate)
	inline Nucleotide RandomNucleotide() const {
		return Nucleotide(RandomBase(), RandomRate(), branchColor, false);
	}

	Tree() : m_nSec(0) {}

	// Evolve the tree
	void Evolve();

	// Add a recombination section to the tree
	bool ProcessTree(NewickNode* pNode);

	template<class itTree>
	bool ProcessTree(itTree itB, itTree itE)
	{
		m_nSec = 0;
		m_map.clear();
		for(itTree it = itB; it!=itE; it++) {
			if(!ProcessTree(*it))
				return false;
		}
		return true;
	}

	const Node::Map& GetMap() const { return m_map; }

	// Align sequences from the tree
	void Align(Alignment &aln, unsigned int uFlags=0);

protected:
	bool ProcessNewickNode(NewickNode* pNode, const std::string &hAnc);
	void Evolve(Sequence &seq, Sequence::const_iterator first, Sequence::const_iterator last, double dTime);
	void Evolve(Node &rNode);
	Tree::Sequence::const_iterator EvolveIndels(Sequence &seq,
		Sequence::const_iterator first, Sequence::const_iterator last,
		double dT);
	size_type NextIndel(double d, double &f);

private:
	int m_nSec;
	Sequence m_vDNASeq;
	Node::Map m_map;
	Node::Map::iterator m_itRoot;
	std::vector<std::string> m_vTips;

	bool m_bRandRootBases, m_bRandRootRates;

	double m_dGamma;
	double m_dIota;

	double m_dTreeScale;

	double m_dOldTime;
	double m_dFreqs[4];
	bitree<double> m_dFreqsCum;
	//Matrix44 m_matSubst;
	Matrix44 m_matTrans;

	bitree<double> m_dTransCum[4];

	Matrix44 m_matV;
	Matrix44 m_matU;
	Matrix44 m_matQ;
	Matrix44 m_matR;
	Vector4  m_vecL;

	std::auto_ptr<IndelModel> m_pInsertionModel;
	std::auto_ptr<IndelModel> m_pDeletionModel;
	double m_dLambdaIns;
	double m_dLambdaDel;
	LinearFunc m_funcRateIns;
	LinearFunc m_funcRateSum;
	int m_uKeepFlank;

	Nucleotide::data_type branchColor;

	residue_factory make_seq;

	std::vector<AlignData> m_vAlnTable;
	typedef std::pair<double, Sequence::size_type> IndelData;
	std::stack<IndelData> m_sInsData, m_sDelData, m_sDelUpData;
};

bool SaveAlignment(std::ostream &rFile, const Tree::Alignment& aln, unsigned int uFlags);

struct AlignData {
	typedef Tree::Sequence Sequence;
	AlignData(const std::string ss) : ssName(ss) {
	}
	AlignData(const AlignData &a) : ssName(a.ssName), seq(a.seq), seqAln(a.seqAln), it(a.it) {
	}
	AlignData & operator=(const AlignData &a) {
		if(this == &a)
			return *this;
		ssName = a.ssName;
		seq = a.seq;
		seqAln = a.seqAln;
		it = a.it;
		return *this;
	}
	std::string ssName;
	const Sequence *seq;
	Sequence seqAln;
	Sequence::const_iterator it;

	struct Printer {
		Printer(const residue_factory &fac) : f(fac) {};

		char operator()(AlignData::Sequence::const_reference r) const {
			static char gaps[] = "-=+";
			return r.is_deleted() ? gaps[r.base()] : f.decode(r.base());
		}
		const residue_factory &f;
	};
};

namespace std {

/*
inline bool operator < (const Tree::Node::Handle & A, const Tree::Node::Handle & B)
{
	return &*A < &*B;
}
*/
}

#endif //DAWG_TREE_H

