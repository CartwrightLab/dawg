#include "node.h"

using namespace std;

////////////////////////////////////////////////////////////
//  class Node
////////////////////////////////////////////////////////////

string Node::ToString() const
{
	string ssTemp = "";
	if(IsClade())
	{
		ssTemp += '(';
		ssTemp += m_pChild->ToString();
		Node* p = m_pChild->m_pSib;
		while(p)
		{
			ssTemp += ',';
			ssTemp += p->ToString();
			p = p->m_pSib;
		}
		ssTemp += ')';
	}
	ssTemp += Label();
	double d = BranchLength();
	if(d != 0.0)
	{
		ssTemp += ':';
		char csBuffer[32];
		sprintf(csBuffer, "%0.10f", BranchLength());
		char *p = &csBuffer[0];
		while (*p) {p++; }
		while(*(--p) == '0' || *p == '.')
			*p = '\0';
		ssTemp += csBuffer;
	}
	return ssTemp;
}

void Node::EvolveSeq()
{
	if(fabs(ScaledLength())< DBL_EPSILON)
		return; // Nothing to evolve
	s_modSubst.Process(this);
	s_modIndel.Process(this);
	//indel formation
	if(g_dLambda < DBL_EPSILON)
		return; // No Indels
	for(double tau = rand_exp(1.0/(g_dLambda*m_seq.size())); tau < dLen;
			tau += rand_exp(1.0/(g_dLambda*m_seq.size())))
		MakeIndel();
}

void Node::Evolve()
{
	if(m_pSib != NULL)
	{
		m_pSib->m_seq = m_seq;
		m_pSib->m_ssGaps = m_ssGaps;
		m_pSib->Evolve();
	}
	EvolveSeq();
	if(m_pChild != NULL)
	{
		m_pChild->m_seq = m_seq;
		m_pChild->m_ssGaps = m_ssGaps;
		m_pChild->Evolve();
	}
}

unsigned int Node::GapPos(unsigned int u) const
{
	string::const_iterator it=m_ssGaps.begin();
	unsigned int v=0;
	//skip leading 'gapspace'
	for(;(*it == '-' || *it == '=') && it != m_ssGaps.end(); ++it)
		v++;
	for(; u && it!=m_ssGaps.end() ;++it)
	{
		v++;
		if(*it != '-' && *it != '=')
			u--;
	}
	return v;
}

void Node::Insert(unsigned int uPos, unsigned int uSize)
{
	//assert(uPos <= m_seq.size());
	m_ssGaps.insert(GapPos(uPos), uSize, '^');

	Seq seq(uSize);
	generate(seq.begin(), seq.end(), Nucleotide::Rand);
	m_seq.insert(m_seq.begin()+uPos, seq.begin(), seq.end());

}
void Node::Delete(unsigned int uPos, unsigned int uSize)
{
	//assert(uPos <= m_seq.size());
	if(uPos+uSize > m_seq.size())
		uSize = m_seq.size()-uPos;
	m_seq.erase(m_seq.begin()+uPos, m_seq.begin()+uPos+uSize);

	uPos = GapPos(uPos);
	while(uSize)
	{
		switch(m_ssGaps[uPos])
		{
		case '=':
		case '-':
			break;
		case '^':
			m_ssGaps[uPos] = '=';
			uSize--;
			break;
		default:
			m_ssGaps[uPos] = '-';
			uSize--;
			break;
		};
		uPos++;
	}
}

////////////////////////////////////////////////////////////
//  class Nucleotide
////////////////////////////////////////////////////////////

bool Nucleotide::Setup(double pFreqs[], double dG, double dI)
{
	s_dGamma = dG;
	s_dIota = dI;
	if(pFreqs[0] < 0 || pFreqs[1] < 0 || pFreqs[2] < 0 || pFreqs[3] < 0)
		return DawgError("Nucleotide frequences need to be positive.");
	memcpy(s_dNucFreqs, pFreqs, 4*sizeof(double));
	s_dNucCumFreqs[0] = pFreqs[0];
	s_dNucCumFreqs[1] = pFreqs[1]+s_dNucCumFreqs[0];
	s_dNucCumFreqs[2] = pFreqs[2]+s_dNucCumFreqs[1];
	s_dNucCumFreqs[3] = 1.0;
	if(s_dNucCumFreqs[0] > 1.0 || s_dNucCumFreqs[1] > 1.0 || s_dNucCumFreqs[2] > 1.0)
		return DawgError("Nucleotide frequences need to sum to 1.");
	return true;
}

Nucleotide Nucleotide::Rand()
{
	Nuc n;
	double d = rand_real();
	if(d <= s_dNucCumFreqs[0])
		n = 0; // A
	else if(d <= s_dNucCumFreqs[1])
		n = 1; // C
	else if(d <= s_dNucCumFreqs[2])
		n = 2; // G
	else
		n = 3; // T

	if(s_dIota > DBL_EPSILON && rand_bool(g_dIota))
		d = 0.0;  // Site Invariant
	else if(s_dGamma > DBL_EPSILON)
		d = rand_gamma1(s_dGamma); // Gamma with mean 1.0 and var of g_dGamma
	else
		d = 1.0;
	return Nucleotide(n, d);
}