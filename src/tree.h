// tree.h - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_TREE_H
#define DAWG_TREE_H

#include "dawg.h"
#include "indel.h"
#include "matrix.h"

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

class Nucleotide
{
public:
	// First two bits specify base
	// Second two bits specify type
	unsigned char    m_ucNuc;
	double m_dRate; // 0.0 means invarant

	Nucleotide() : m_ucNuc(0xF), m_dRate(1.0) { }
	Nucleotide(unsigned char nuc, double rate) : m_ucNuc(nuc), m_dRate(rate) { }

	static const unsigned char MaskBase		= 0x3; // 0011
	static const unsigned char MaskType		= 0xC; // 1100
	static const unsigned char MaskDel		= 0x8; // 1000
	static const unsigned char MaskIns		= 0x4; // 0100
	static const unsigned char TypeRoot		= 0x0; // 0000
	static const unsigned char TypeIns		= 0x4; // 0100
	static const unsigned char TypeDel		= 0x8; // 1000
	static const unsigned char TypeDelIns	= 0xC; // 1100
	
	inline unsigned char GetBase() const { return m_ucNuc & MaskBase; }
	inline unsigned char GetType() const { return m_ucNuc & MaskType; }
	inline void SetBase(unsigned char uc) { m_ucNuc =  (uc & MaskBase) | (m_ucNuc & MaskType); }
	inline void SetType(unsigned char uc) { m_ucNuc =  (uc & MaskType) | (m_ucNuc & MaskBase); }
	inline void SetNuc(unsigned char ucB, unsigned char ucT)
		{ m_ucNuc =  (ucB & MaskBase) | (ucT & MaskType); }
	inline bool IsType(unsigned char uc) const { return (GetType() == uc); }
	inline bool IsDeletion() const { return ((m_ucNuc & MaskDel) == MaskDel); }
	inline bool IsInsertion() const { return ((m_ucNuc & MaskIns) == MaskIns); }

	inline bool FromChar(char ch)
	{
		switch(ch)
		{
		case 'A':
		case 'a':
			m_ucNuc = 0;
			return true;
		case 'C':
		case 'c':
			m_ucNuc = 1;
			return true;
		case 'G':
		case 'g':
			m_ucNuc = 2;
			return true;
		case 'T':
		case 't':
			m_ucNuc = 3;
			return true;
		default:
			return false;
		}
	}
	inline char ToChar() const
	{
		char csNuc[]	= "ACGT";
		char csType[]	= " +-=";
		return IsType(TypeRoot) ? csNuc[GetBase()] : csType[GetType() >> 2];
	}
};

class Sequence : public std::vector<Nucleotide>
{
public:
	typedef std::vector<Nucleotide> Base;
	Sequence() : m_uLength(0) { }
	Sequence(unsigned long uSize) : Base(uSize, Nucleotide(0xF, -1.0))
	{
		m_uLength = uSize;
	}
	unsigned long SeqLength() const { return m_uLength; }

	const_iterator SeqPos(unsigned long uPos) const;
	iterator SeqPos(unsigned long uPos);

	unsigned long Insertion(iterator itPos, const_iterator itBegin, const_iterator itEnd);
	unsigned long Deletion(iterator itBegin, unsigned long uSize);

	void Append(const Sequence &seq);

	void ToString(std::string &ss) const;

private:
	unsigned long m_uLength;
};

class Tree
{
public:
	class Node
	{
	public:
		typedef std::map<std::string, Tree::Node> Map;
		typedef Map::iterator Handle;
		std::vector<Sequence> m_vSections;
		std::vector<Handle> m_vAncestors;
		std::map<Handle, double> m_mBranchLens;
		bool m_bTouched;
		Node() : m_bTouched(false) { }
		void Flatten(Sequence& seq) const;
		unsigned long SeqLength() const;

		typedef std::pair<std::vector<Sequence>::iterator, Sequence::iterator> iterator;
		typedef std::pair<std::vector<Sequence>::const_iterator, Sequence::const_iterator> const_iterator;

		iterator SeqPos(unsigned long uPos);
		const_iterator SeqPos(unsigned long uPos) const;
	};

	typedef std::map<std::string, std::string> Alignment;
	
	bool SetupEvolution(double pFreqs[], double pSubs[],
		const IndelModel::Params& rIns, const IndelModel::Params& rDel,
		unsigned long uWidth, const std::vector<double> &vdGamma,
		const std::vector<double> &vdIota, const std::vector<double> &vdScale, double dTreeScale);
	bool SetupRoot(const std::vector<std::string> &vSeqs, const std::vector<unsigned long> &vData,
		const std::vector<std::vector<double> > &vRates);

	double RandomRate(unsigned long uPos) const;
	unsigned char RandomNuc() const;
	Nucleotide RandomNucleotide(unsigned long uPos) const
		{ return Nucleotide(RandomNuc(), RandomRate(uPos)); }

	Tree() : m_nSec(0), m_uWidth(1) {}
	
	inline unsigned long BlockTrim(unsigned long u) { return u - u%m_uWidth; }

	void Evolve();
	void ProcessTree(NewickNode* pNode);
	
	const Node::Map& GetMap() const { return m_map; }

	void Align(Alignment &aln, bool bGapPlus, bool bGapSingleChar) const;

protected:
	void ProcessNewickNode(NewickNode* pNode, Node::Handle hAnc);
	void Evolve(Node &rNode, double dTime);
	void Evolve(Node &rNode);

private:
	int m_nSec;
	std::vector< Sequence > m_vDNASeqs;
	Node::Map m_map;

	unsigned long m_uWidth;
	std::vector<double> m_vdScale;
	std::vector<double> m_vdGamma;
	std::vector<double> m_vdIota;

	double m_dTreeScale;

	double m_dOldTime;
	double m_dFreqs[4];
	double m_dNucCumFreqs[4];
	Matrix44 m_matSubst;
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
};

inline bool operator < (const Tree::Node::Handle & A, const Tree::Node::Handle & B)
{
	return &*A < &*B;
}

bool SaveAlignment(std::ostream &rFile, const Tree::Alignment& aln);

#endif //DAWG_TREE_H

