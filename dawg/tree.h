#ifndef DAWG_TREE_H
#define DAWG_TREE_H

#include "dawg.h"
#include "indel.h"
#include "matrix.h"

class Nucleotide;

class Sequence
{
public:
	Sequence();
	Sequence(unsigned long uSize);
	Sequence(std::string ssDNA);

	unsigned long GapPos(unsigned long uPos) const;
	
	bool Insert(unsigned long uPos, unsigned long uSize);
	bool Delete(unsigned long uPos, unsigned long uSize);

	// Access the dna sequence
	Nucleotide& operator [](unsigned long uPos) { return m_vDNA[uPos]; }
	const Nucleotide& operator [](unsigned long uPos) const { return m_vDNA[uPos]; }

	unsigned long Length() const { return m_vDNA.size(); }

private:
	std::vector<Nucleotide> m_vDNA;
	std::vector<char> m_vHistory;
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
		bool Insert(unsigned long uPos, unsigned long uSize);
		bool Delete(unsigned long uPos, unsigned long uSize);
	};
	
	bool SetupSubst(double pFreqs[], double pSubs[]);
	bool SetupIndel(const IndelModel::Params& rIns, const IndelModel::Params& rDel);
	
	Tree() : m_nSec(0) {}
	
	void Evolve();
	void Process(NewickNode* pNode);
	
	const Node::Map& GetMap() const { return m_map; }

protected:
	void ProcessNode(NewickNode* pNode);
	//void Touch(node& rNode);
	void Evolve(Node &rNode, double dTime);

private:
	int m_nSec;
	Node::Map m_map;
	std::vector<unsigned int> m_vuSecLength;
	double m_dScale;

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

class Nucleotide
{
public:
	typedef unsigned char Nuc;

	Nucleotide() : m_nuc(5), m_dRate(1.0) { }
	Nucleotide(Nuc nuc, double rate) : m_nuc(nuc), m_dRate(rate) { }

	double m_dRate; // 0.0 means invarant
	Nuc    m_nuc;

	static Nucleotide Rand();
	static bool Setup(double pFreqs[], double dG, double dI);
	
	static double Iota();
	static double Gamma();

protected:
	static double s_dNucCumFreqs[4];
	static double s_dNucFreqs[4];
	static double s_dGamma;
	static double s_dIota;
};

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

bool SaveSequences(std::ostream &rFile, const Node *arTrees[], unsigned int uSize);

#endif //DAWG_TREE_H