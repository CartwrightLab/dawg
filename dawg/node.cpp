#include "dawg.h"
#include "node.h"
#include "rand.h"

using namespace std;

unsigned long hash_adler32(const char* cs);

////////////////////////////////////////////////////////////
//  class Node
////////////////////////////////////////////////////////////
Node::Node(Node* p, const char *cs, double d) : m_dLen(d), m_pSub(p)
{
	m_csLabel[1023] = '\0';
	if(cs)
		strncpy(m_csLabel, cs, 1023);
	else
		MakeName();
}

void Node::MakeName()
{
	std::vector<std::string> v;
	Node* p = m_pSub.get();
	while(p != NULL)
	{
		v.push_back(p->m_csLabel);
		p = p->m_pSib.get();
	}
	std::sort(v.begin(), v.end());
	std::string ss;
	for(std::vector<std::string>::iterator it = v.begin();
		it != v.end(); ++it)
	{
		ss += *it;
		ss += ",";
	}
	m_csLabel[0] = '!'; 
	_ltoa(hash_adler32(ss.c_str()), m_csLabel+1, 16);
}

////////////////////////////////////////////////////////////
//  class Tree
////////////////////////////////////////////////////////////

void Tree::Process(Node* pNode)
{
	ProcessNode(pNode);
	m_nSec++;
}

void Tree::ProcessNode(Node* pNode)
{
	static vector<Map::iterator> vStack;

	if(pNode->m_pSib.get())
		ProcessNode(pNode->m_pSib.get());

	Section sec;
	sec.dBranchLen = pNode->m_dLen;
	if(vStack.empty())
		sec.itAncestor = m_map.end();
	else
		sec.itAncestor = vStack.back();
	node& rn = m_map[pNode->m_csLabel];
	for(node::const_iterator it = rn.begin(); it != rn.end(); ++it)
	{
		if(sec.itAncestor == it->itAncestor)
			sec.dBranchLen = it->dBranchLen;
	}
	rn.resize(m_nSec+1, Section(0.0, m_map.end()));
	rn[m_nSec] = sec;
			
	if(pNode->m_pSub.get())
	{
		vStack.push_back(m_map.find(pNode->m_csLabel));
		ProcessNode(pNode->m_pSub.get());
		vStack.pop_back();
	}
}

////////////////////////////////////////////////////////////
//  class Node
////////////////////////////////////////////////////////////

//double Node::s_dScale;
//IndelProcessor Node::s_procIndel;
//SubstProcessor Node::s_procSubst;
//
//string Node::ToString() const
//{
//	string ssTemp = "";
//	if(IsClade())
//	{
//		ssTemp += '(';
//		ssTemp += m_pChild->ToString();
//		Node* p = m_pChild->m_pSib.get();
//		while(p)
//		{
//			ssTemp += ',';
//			ssTemp += p->ToString();
//			p = p->m_pSib.get();
//		}
//		ssTemp += ')';
//	}
//	ssTemp += Label();
//	double d = BranchLength();
//	if(d != 0.0)
//	{
//		ssTemp += ':';
//		char csBuffer[32];
//		sprintf(csBuffer, "%0.10f", BranchLength());
//		char *p = &csBuffer[0];
//		while (*p) {p++; }
//		while(*(--p) == '0' || *p == '.')
//			*p = '\0';
//		ssTemp += csBuffer;
//	}
//	return ssTemp;
//}
//
//void Node::EvolveSeq()
//{
//	if(fabs(ScaledLength())< DBL_EPSILON)
//		return; // Nothing to evolve
//	s_procSubst.Process(this);
//	s_procIndel.Process(this);
//}
//
//void Node::Evolve()
//{
//	if(m_pSib.get() != NULL)
//	{
//		m_pSib->m_seq = m_seq;
//		m_pSib->m_ssGaps = m_ssGaps;
//		m_pSib->Evolve();
//	}
//	EvolveSeq();
//	if(m_pChild.get() != NULL)
//	{
//		m_pChild->m_seq = m_seq;
//		m_pChild->m_ssGaps = m_ssGaps;
//		m_pChild->Evolve();
//	}
//}
//
//unsigned int Node::GapPos(unsigned int u) const
//{
//	string::const_iterator it=m_ssGaps.begin();
//	unsigned int v=0;
//	//skip leading 'gapspace'
//	for(;(*it == '-' || *it == '=') && it != m_ssGaps.end(); ++it)
//		v++;
//	for(; u && it!=m_ssGaps.end() ;++it)
//	{
//		v++;
//		if(*it != '-' && *it != '=')
//			u--;
//	}
//	return v;
//}
//
//void Node::Insert(unsigned int uPos, unsigned int uSize)
//{
//	//assert(uPos <= m_seq.size());
//	m_ssGaps.insert(GapPos(uPos), uSize, '^');
//
//	Seq seq(uSize);
//	generate(seq.begin(), seq.end(), Nucleotide::Rand);
//	m_seq.insert(m_seq.begin()+uPos, seq.begin(), seq.end());
//
//}
//void Node::Delete(unsigned int uPos, unsigned int uSize)
//{
//	//assert(uPos <= m_seq.size());
//	uPos++;
//	unsigned int uStart = (uSize < uPos)? uPos-uSize : 0u;
//	unsigned int uEnd = (m_seq.size() > uPos) ? uPos : m_seq.size();
//
//	m_seq.erase(m_seq.begin()+uStart, m_seq.begin()+uEnd);
//
//	uSize = uEnd-uStart;
//	uPos = GapPos(uStart);
//	while(uSize)
//	{
//		switch(m_ssGaps[uPos])
//		{
//		case '=':
//		case '-':
//			break;
//		case '^':
//			m_ssGaps[uPos] = '=';
//			uSize--;
//			break;
//		default:
//			m_ssGaps[uPos] = '-';
//			uSize--;
//			break;
//		};
//		uPos++;
//	}
//}

////////////////////////////////////////////////////////////
//  class Nucleotide
////////////////////////////////////////////////////////////

double Nucleotide::s_dNucCumFreqs[4] = {0.25, 0.50, 0.75, 1.0};
double Nucleotide::s_dNucFreqs[4] = {0.25, 0.25, 0.25, 0.25};
double Nucleotide::s_dGamma = 0.0;
double Nucleotide::s_dIota = 0.0;

double Nucleotide::Gamma()
{
	return Nucleotide::s_dGamma;
}

double Nucleotide::Iota()
{
	return Nucleotide::s_dIota;
}

bool Nucleotide::Setup(double pFreqs[], double dG, double dI)
{
	s_dGamma = dG;
	s_dIota = dI;
	if(s_dGamma < 0.0)
		return DawgError("Invalid Gamma, \"%f\".  Gamma must be positive.", s_dGamma);
	else if(0.0 > s_dIota || s_dIota > 1.0)
		return DawgError("Invalid Iota, \"%f\".  Iota must be a probability.", s_dIota);

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

	if(s_dIota > DBL_EPSILON && rand_bool(s_dIota))
		d = 0.0;  // Site Invariant
	else if(s_dGamma > DBL_EPSILON)
		d = rand_gamma1(s_dGamma); // Gamma with mean 1.0 and var of g_dGamma
	else
		d = 1.0;
	return Nucleotide(n, d);
}

////////////////////////////////////////////////////////////
//  Utility Functions
////////////////////////////////////////////////////////////

unsigned long hash_adler32(const char* cs)
{	// From the Crypt++ Library
	const unsigned long BASE = 65521;
	unsigned long s1 = 1;
	unsigned long s2 = 0;
	size_t len = strlen(cs);

	if (len % 8 != 0)
	{
		do
		{
			s1 += *cs++;
			s2 += s1;
			len--;
		} while (len % 8 != 0);

		if (s1 >= BASE)
			s1 -= BASE;
		s2 %= BASE;
	}
	while (len > 0)
	{
		s1 += cs[0]; s2 += s1;
		s1 += cs[1]; s2 += s1;
		s1 += cs[2]; s2 += s1;
		s1 += cs[3]; s2 += s1;
		s1 += cs[4]; s2 += s1;
		s1 += cs[5]; s2 += s1;
		s1 += cs[6]; s2 += s1;
		s1 += cs[7]; s2 += s1;

		len -= 8;
		cs += 8;

		if (s1 >= BASE)
			s1 -= BASE;
		if (len % 0x8000 == 0)
			s2 %= BASE;
	}
	return (s1 << 16) | (s2 & 0xFFFF);
}