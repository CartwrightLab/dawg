#ifndef DAWG_NODE_H
#define DAWG_NODE_H

#include <string>
#include <vector>
#include <memory>

class Nucleotide;
typedef std::vector<Nucleotide> Seq;

class Node
{
public:
	Node(const char* csName = "", Node * pChild = NULL) : m_ssLabel(csName),
		m_pChild(pChild), m_pSib(NULL), m_dBranchLength(0.0) { }

	bool IsTaxon() const { return m_pChild == NULL; }
	bool IsClade() const { return m_pChild != NULL; }

	void AddChild(Node* pNode)
	{
		pNode->m_pSib = m_pChild;
		m_pChild = pNode;
	}
	bool IsChild(Node* pNode)
	{
		Node *p = m_pChild;
		while(p != pNode && p)
			p = p->m_pSib;
		return p != NULL;
	}
	void AddSib(Node* pNode)
	{
		pNode->m_pSib = m_pSib;
		m_pSib = pNode;
	}
	
	const std::string& Label() const { return m_ssLabel; }
	void BranchLength(double d) { m_dBranchLength = d; }
	double BranchLength() const { return m_dBranchLength; }
	double ScaledLength() const { return s_dScale*m_dBranchLength; }

	std::string ToString() const;

	const Seq& Sequence() const { return m_seq; }
	Seq& Sequence() { return m_seq; }

	const Node *Child() const {return m_pChild; }
	const Node *Sibling() const {return m_pSib; }

	void Evolve();

	const std::string& Gaps() const { return m_ssGaps; }
	void  ResetGaps() { m_ssGaps.assign(m_seq.size(), '%'); }
	unsigned int GapPos(unsigned int uPos) const;
	void Insert(unsigned int uPos, unsigned int uSize);
	void Delete(unsigned int uPos, unsigned int uSize);
	
	static void Scale(double d) { s_dScale = d; }
	static double Scale() { return s_dScale; }
	
	static IndelProcessor s_procIndel;
	static SubstProcessor s_procSubst

protected:
	std::auto_ptr<Node> m_pSib;
	std::auto_ptr<Node> m_pChild;
	double m_dBranchLength;
	std::string m_ssLabel;

	Seq m_seq;
	std::string m_ssGaps;

private:
	Node(const Node&);
	Node& operator=(const Node&);

	void EvolveSeq();

	static double s_dScale;
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

protected:
	static s_dNucCumFreqs[4];
	static s_dNucFreqs[4];
	static s_dGamma;
	static s_dIota;
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

#endif //DAWG_NODE_H