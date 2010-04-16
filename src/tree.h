// tree.h - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_TREE_H
#define DAWG_TREE_H

#include "dawg.h"
#include "indel.h"
#include "matrix.h"

#include <dawg/residue.h>

#include "bitree.h"
#include "rand.h"

using namespace dawg;

#include <list>
#include <stack>
#include <vector>

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

struct AlignData;

// The recombinant tree data structure
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
		if(m_bIota && rand_bool(m_dIota))
			return Nucleotide::rate_type(0.0);  // Site Invariant
		return static_cast<Nucleotide::rate_type>(rand_gamma1(m_dGamma)); // Gamma with mean 1.0 and var of m_dGamma
	}
	bool AreRatesConstant() const {
		return m_bConstRates;
	}
	// Draw a random base from the evolutionary parameters
	inline Nucleotide::data_type RandomBase() const {
		return m_dFreqsCum(rand_real());
	}

	// Draw a random nucleotide (base and rate)
	inline Nucleotide RandomNucleotide() const {
		return Nucleotide(RandomBase(), RandomRate(), branchColor);
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

	bool m_bRandRootBases, m_bRandRootRates, m_bConstRates, m_bIota;

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
	AlignData(const Sequence *xseq, Tree::Alignment::mapped_type *xstr) :
		seq(xseq), str(xstr), it(xseq->begin()) {
	}
	const Sequence *seq;
	Tree::Alignment::mapped_type *str;
	Sequence::const_iterator it;
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

