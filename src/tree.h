#ifndef DAWG_TREE_H
#define DAWG_TREE_H

#include "dawg.h"
#include "indel.h"
#include "matrix.h"

class Nucleotide
{
public:
	typedef unsigned char Nuc;

	Nucleotide() : m_nuc(5), m_dRate(1.0) { }
	Nucleotide(Nuc nuc, double rate) : m_nuc(nuc), m_dRate(rate) { }

	double m_dRate; // 0.0 means invarant
	Nuc    m_nuc;
};

class Sequence
{
public:
	typedef std::vector<Nucleotide> DNAVec;
	typedef std::vector<char> HistoryVec;

	Sequence();
	Sequence(const DNAVec &dna);
	unsigned long GapPos(unsigned long uPos) const;
	
	unsigned long Insert(unsigned long uPos, DNAVec::const_iterator itBegin, DNAVec::const_iterator itEnd);
	unsigned long Delete(unsigned long uPos, unsigned long uSize);

	// Access the dna sequence
	Nucleotide& operator [](unsigned long uPos) { return m_vDNA[uPos]; }
	const Nucleotide& operator [](unsigned long uPos) const { return m_vDNA[uPos]; }

	unsigned long Length() const { return m_vDNA.size(); }
	
	const DNAVec& DNA() const { return m_vDNA; }
	const HistoryVec& History() const {return m_vHistory; } 
	
	void Append(const Sequence& seq);

	void ResetHistory();

private:
	DNAVec		m_vDNA;
	HistoryVec	m_vHistory;
};

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
		unsigned long SeqLength() const;
		unsigned long Insert(unsigned long uPos,
			Sequence::DNAVec::const_iterator itBegin,
			Sequence::DNAVec::const_iterator itEnd);
		unsigned long Delete(unsigned long uPos, unsigned long uSize);
		void Flatten(Sequence& seq) const;
	};
	typedef std::map<std::string, std::string> Alignment;
	
	bool SetupEvolution(double pFreqs[], double pSubs[],
		const IndelModel::Params& rIns, const IndelModel::Params& rDel,
		double dGamma, double dIota, double dScale);
	bool SetupRoot(const std::vector<std::string> &vSeqs, const std::vector<int> &vData,
		const std::vector<std::vector<double> > &vRates);

	double RandomRate() const;
	Nucleotide::Nuc RandomNuc() const;
	Nucleotide RandomNucleotide() const { return Nucleotide(RandomNuc(), RandomRate()); }

	Tree() : m_nSec(0) {}
	
	void Evolve();
	void ProcessTree(NewickNode* pNode);
	
	const Node::Map& GetMap() const { return m_map; }

	void Align(Alignment &aln) const;

protected:
	void ProcessNewickNode(NewickNode* pNode, Node::Handle hAnc);
	void Evolve(Node &rNode, double dTime);
	void Evolve(Node &rNode);

private:
	int m_nSec;
	std::vector< Sequence::DNAVec > m_vDNASeqs;
	Node::Map m_map;
	double m_dScale;

	double m_dGamma;
	double m_dIota;

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

inline Nucleotide::Nuc CharToNuc(char ch)
{
	switch(ch)
	{
		case 'a':
		case 'A': return 0;
		case 'c':
		case 'C': return 1;
		case 'g':
		case 'G': return 2;
		case 't':
		case 'T': return 3;
		case '-': return 4;
		case '?': return 5;
		default : return (Nucleotide::Nuc)-1;
	}
}

inline char NucToChar(Nucleotide::Nuc n)
{
	static char cs[] = "ACGT-?";
	return (n > 4) ? '?' : cs[n];
}

bool SaveAlignment(std::ostream &rFile, const Tree::Alignment& aln);

#endif //DAWG_TREE_H

