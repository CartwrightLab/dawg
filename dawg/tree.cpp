#include "dawg.h"
#include "tree.h"
#include "rand.h"

using namespace std;

////////////////////////////////////////////////////////////
//  class NewickNode
////////////////////////////////////////////////////////////
NewickNode::NewickNode(NewickNode* p, const char *cs, double d) : m_dLen(d), m_pSub(p)
{
	if(cs)
		m_ssLabel = cs;
	else
		MakeName();
}

void NewickNode::MakeName()
{
	vector<const std::string*> v;
	NewickNode* p = m_pSub.get();
	while(p != NULL)
	{
		v.push_back(&p->m_ssLabel);
		p = p->m_pSib.get();
	}
	std::sort(v.begin(), v.end());
	m_ssLabel = "(";
	m_ssLabel += *v.front();
	for(vector<const std::string*>::iterator it = v.begin()+1; it != v.end(); ++it)
	{
		m_ssLabel += ",";
		m_ssLabel += **it;	
	}
	m_ssLabel += ")";
}

////////////////////////////////////////////////////////////
//  class Sequence
////////////////////////////////////////////////////////////
Sequence::Sequence()
{

}

Sequence::Sequence(unsigned long uSize) : m_vDNA(uSize), m_vHistory(uSize, '.')
{
	generate(m_vDNA.begin(), m_vDNA.end(), Nucleotide::Rand);
}

Sequence::Sequence(std::string ssDNA) : m_vDNA(ssDNA.length()), m_vHistory(ssDNA.length(), '.')
{
	generate(m_vDNA.begin(), m_vDNA.end(), Nucleotide::Rand);
	for(unsigned long u = 0; u < m_vDNA.size(); ++u)
		m_vDNA[u].m_nuc = CharToNuc(ssDNA[u]);
}

unsigned long Sequence::GapPos(unsigned long uPos) const
{
	vector<char>::const_iterator it=m_vHistory.begin();
	unsigned long v=0;
	//skip leading 'gapspace'
	for(;(*it == '-' || *it == '=') && it != m_vHistory.end(); ++it)
		v++;
	for(; uPos && it != m_vHistory.end(); ++it)
	{
		v++;
		if(*it != '-' && *it != '=')
			uPos--;
	}
	return v;
}

unsigned long Sequence::Insert(unsigned long uPos, unsigned long uSize)
{
	if(uSize == 0 || uPos > m_vDNA.size())
		return 0;
	m_vHistory.insert(m_vHistory.begin()+GapPos(uPos), uSize, '+');
	vector<Nucleotide> seq(uSize);
	generate(seq.begin(), seq.end(), Nucleotide::Rand);
	m_vDNA.insert(m_vDNA.begin()+uPos, seq.begin(), seq.end());
	return uSize;
}

unsigned long Sequence::Delete(unsigned long uPos, unsigned long uSize)
{
	if(uSize == 0)
		return 0;
	uPos++;
	unsigned long uStart = (uSize < uPos)? uPos-uSize : 0u;
	unsigned long uEnd = (m_vDNA.size() > uPos) ? uPos : m_vDNA.size();
	
	m_vDNA.erase(m_vDNA.begin()+uStart, m_vDNA.begin()+uEnd);

	unsigned long uTemp = uSize = uEnd-uStart;
	uPos = GapPos(uStart);
	while(uTemp)
	{
		switch(m_vHistory[uPos])
		{
		case '=':
		case '-':
			break;
		case '+':
			m_vHistory[uPos] = '=';
			uTemp--;
			break;
		default:
			m_vHistory[uPos] = '-';
			uTemp--;
			break;
		};
		uPos++;
	}
	return uSize;
}


////////////////////////////////////////////////////////////
//  class Tree::Node
////////////////////////////////////////////////////////////

unsigned long Tree::Node::SeqLength() const
{
	unsigned long uRet = 0;
	for(vector<Sequence>::const_iterator it = m_vSections.begin(); it != m_vSections.end(); ++it)
		uRet += it->Length();
	return uRet;
}

unsigned long Tree::Node::Insert(unsigned long uPos, unsigned long uSize)
{
	if(uSize == 0)
		return 0;	
	vector<Sequence>::iterator it;
	for(it = m_vSections.begin(); it != m_vSections.end() && uPos > it->Length(); ++it)
			uPos -= it->Length();
	if( it == m_vSections.end())
		return 0;
	else if(uPos < it->Length())
		return it->Insert(uPos, uSize);
	else
	{
		// Insertion occurs at a gap between sections
		// Randomly allocate uSize among the sections
		vector<unsigned long> vTemp;
		vTemp.push_back(0);
		vTemp.push_back(rand_ulong(uSize));
		for(vector<Sequence>::iterator jt = it; jt != m_vSections.end() && jt->Length(); ++jt)
			vTemp.push_back(rand_ulong(uSize));
		vTemp.push_back(uSize);
		sort(vTemp.begin(), vTemp.end());
		unsigned long uTemp = it->Insert(uPos, vTemp[1]-vTemp[0]);
		int i=1;
		for(vector<Sequence>::iterator jt = it+1; jt != m_vSections.end() && jt->Length(); ++jt, ++i)
			uTemp += jt->Insert(0, vTemp[i+1]-vTemp[i]);
		return uTemp;
	}
}

unsigned long Tree::Node::Delete(unsigned long uPos, unsigned long uSize)
{
	if(uSize == 0)
		return 0;
	vector<Sequence>::iterator it;
	for(it = m_vSections.begin(); it != m_vSections.end() && uPos > it->Length(); ++it)
			uPos -= it->Length();
	if( it == m_vSections.end())
		return 0;
	return it->Delete(uPos, uSize);
}

////////////////////////////////////////////////////////////
//  class Tree
////////////////////////////////////////////////////////////

void Tree::ProcessTree(NewickNode* pNode)
{
	ProcessNewickNode(pNode);
	m_nSec++;
}

void Tree::ProcessNewickNode(NewickNode* pNode)
{
	// NOT THREAD SAFE!!!!!!
	static vector<Node::Handle> vStack;

	if(pNode->m_pSib.get())
		ProcessNewickNode(pNode->m_pSib.get());
	
	Node::Handle hAnc = (vStack.size()) ? vStack.back() : m_map.end();
	Node& node = m_map[pNode->m_ssLabel];
	node.m_mBranchLens[hAnc] = pNode->m_dLen;
	node.m_vAncestors.resize(m_nSec+1, m_map.end());
	node.m_vAncestors[m_nSec] = hAnc;
			
	if(pNode->m_pSub.get())
	{
		vStack.push_back(m_map.find(pNode->m_ssLabel));
		ProcessNewickNode(pNode->m_pSub.get());
		vStack.pop_back();
	}
}

void Tree::Evolve()
{
	for(Node::Handle it = m_map.begin(); it != m_map.end(); ++it)
		Evolve(it->second);
}

void Tree::Evolve(Node &rNode)
{
	if(rNode.m_bTouched)
		return;
	rNode.m_bTouched = true;
	map<Node::Handle, Node> mapSeqs;
	for(unsigned long a = 0; a < rNode.m_vAncestors.size(); ++a)
	{
		if(rNode.m_vAncestors[a] == m_map.end())
		{
			// no ancestor need to create one
		}
		else
		{
			if(mapSeqs.find(rNode.m_vAncestors[a]) == mapSeqs.end())
			{
				Evolve(rNode.m_vAncestors[a]->second);
				mapSeqs[rNode.m_vAncestors[a]] = rNode.m_vAncestors[a]->second;
				Evolve(mapSeqs[rNode.m_vAncestors[a]], rNode.m_mBranchLens[rNode.m_vAncestors[a]]);
			}
			rNode.m_vSections.push_back(mapSeqs[rNode.m_vAncestors[a]].m_vSections[a]);
		}
	}
}

void Tree::Evolve(Node &rNode, double dTime)
{
	dTime = fabs(dTime);
	if(dTime < DBL_EPSILON)
		return; // Nothing to evolve
	
	// Substitutions
	for(vector<Sequence>::iterator it = rNode.m_vSections.begin(); it != rNode.m_vSections.end(); ++it)
		for(unsigned long u = 0; u < it->Length(); ++it)
		{
			if((*it)[u].m_dRate < DBL_EPSILON)
				continue; // Invariant Site
			double dTemp = dTime*(*it)[u].m_dRate;
			if(dTemp != m_dOldTime)
			{
				m_dOldTime = dTemp;
				Vector4  vec;
				vec[0] = exp(dTemp*m_vecL[0]);
				vec[1] = exp(dTemp*m_vecL[1]);
				vec[2] = exp(dTemp*m_vecL[2]);
				vec[3] = exp(dTemp*m_vecL[3]);
				Matrix44 mat; mat.Scale(vec, m_matU);
				m_matSubst.Multiply(m_matV, mat);
				for(Matrix44::Pos i=0;i<4;++i)
				{
					m_matSubst(i,1) += m_matSubst(i,0);
					m_matSubst(i,2) += m_matSubst(i,1);
					//m_matSubst(i,3) = 1.0;
				}
			}
			dTemp = rand_real();
			if(dTemp <= m_matSubst((*it)[u].m_nuc, 0))
				(*it)[u].m_nuc = 0;
			else if(dTemp <= m_matSubst((*it)[u].m_nuc, 1))
				(*it)[u].m_nuc = 1;
			else if(dTemp <= m_matSubst((*it)[u].m_nuc, 2))
				(*it)[u].m_nuc = 2;
			else
				(*it)[u].m_nuc = 3;
		}
	// Indels
	if(m_funcRateSum(0.0) < DBL_EPSILON)
		return;

	unsigned long uLength = rNode.SeqLength();
	double dLength = (double)uLength;
	double dW = 1.0/m_funcRateSum(dLength);
	double dt = rand_exp(dW);
	while(dt <= dTime)
	{
		if(rand_bool(m_funcRateIns(dLength)*dW))
		{
			//Insertion
			uLength += rNode.Insert(rand_ulong(uLength), m_pInsertionModel->RandSize());
		}
		else
		{
			//Deletion
			unsigned long ul = m_pDeletionModel->RandSize();
			uLength -= rNode.Delete(rand_ulong(uLength+ul-1), ul);
		}
		dLength = (double)uLength;
		dW = 1.0/m_funcRateSum(dLength);
		dt += rand_exp(dW);
	}
}

bool Tree::SetupSubst(double pFreqs[], double pSubs[])
{
	if(pFreqs[0] < 0.0 || pFreqs[1] < 0.0 || pFreqs[2] < 0.0 || pFreqs[3] < 0.0)
		return DawgError("Nucleotide frequences need to be positive.");
	pFreqs[3] = 1.0-pFreqs[0]-pFreqs[1]-pFreqs[2];
	if( pFreqs[3] < 0.0 )
		return DawgError("Nucleotide frequencies need to sum to 1.0.");
	if(pSubs[0] < 0.0 || pSubs[1] < 0.0 || pSubs[2] < 0.0
		|| pSubs[3] < 0.0 || pSubs[4] < 0.0 || pSubs[5] < 0.0)
		return DawgError("Substitution rates need to be positive.");

	// Setup Symetric Matrix
	Matrix44 matQ(Matrix44::s_Zero);
	matQ(0,1) = matQ(1,0) = pSubs[0]; //A-C
	matQ(0,2) = matQ(2,0) = pSubs[1]; //A-G
	matQ(0,3) = matQ(3,0) = pSubs[2]; //A-T
	matQ(1,2) = matQ(2,1) = pSubs[3]; //C-G
	matQ(1,3) = matQ(3,1) = pSubs[4]; //C-T
	matQ(2,3) = matQ(3,2) = pSubs[5]; //G-T
	
	// Store Rate Matrix
	m_matR = matQ;
	
	// Create GTR Genetating Matrix
	Vector4 vecF(pFreqs);	
	matQ.Scale(matQ, vecF);

	// Scale such that the total rate of substitution is equal to one
	matQ(0,0) = -(matQ(0,1)+matQ(0,2)+matQ(0,3));
	matQ(1,1) = -(matQ(1,0)+matQ(1,2)+matQ(1,3));
	matQ(2,2) = -(matQ(2,0)+matQ(2,1)+matQ(2,3));
	matQ(3,3) = -(matQ(3,0)+matQ(3,1)+matQ(3,2));
	matQ.Scale(matQ, -1.0/((1.0-Nucleotide::Iota())*(vecF[0]*matQ(0,0)+vecF[1]*matQ(1,1)+
		vecF[2]*matQ(2,2)+vecF[3]*matQ(3,3))));
	
	// Store Scaled Q Matrix
	m_matQ = matQ;

	// Make Q a symetric matrix again
	Vector4 vecD, vecE;  //D*E=I
	for(Matrix44::Pos i=0;i<4;++i)
	{
		vecD[i] = sqrt(vecF[i]);
		vecE[i] = 1.0/vecD[i];
	}
	matQ.Scale(matQ, vecE);
	matQ.Scale(vecD, matQ);

	//Find EigenSystem using Jacobian Transformations
	int nRet = EigenSystem(matQ, m_vecL, m_matV);
	if(nRet == -1)
		return DawgError("Eigensystem failed to converge.");
	m_matU = m_matV;
	m_matV.Scale(vecE, m_matV);
	m_matU.Transpose();
	m_matU.Scale(m_matU, vecD);
	return true;
}

bool Tree::SetupIndel(const IndelModel::Params& rIns, const IndelModel::Params& rDel)
{
	m_dLambdaIns = rIns.dLambda;
	if(m_dLambdaIns < 0.0)
		return DawgError("Lambda (Ins) must not be negative.");
	if(rIns.ssModel == "NB")
	{
		try {m_pInsertionModel.reset(new NegBnModel(rIns.vdModel));}
			catch(...) {return DawgError("Insertion model parameters not specified correctly.");}
	}
	else
	{
		try {m_pInsertionModel.reset(new UserModel(rIns.vdModel));}
			catch(...) {return DawgError("Insertion model parameters not specified correctly.");}
	}

	m_dLambdaDel = rDel.dLambda;
	if(m_dLambdaDel < 0.0)
		return DawgError("Lambda (Del) must not be negative.");
	if(rDel.ssModel == "NB")
	{
		try {m_pDeletionModel.reset(new NegBnModel(rDel.vdModel));}
			catch(...) {return DawgError("Deletion model parameters not specified correctly.");}
	}
	else
	{
		try {m_pDeletionModel.reset(new UserModel(rDel.vdModel));}
			catch(...) {return DawgError("Deletion model parameters not specified correctly.");}
	}    
	m_funcRateIns.m = m_dLambdaIns;
	m_funcRateIns.b = m_dLambdaIns;
	m_funcRateSum.m = m_dLambdaDel+m_dLambdaIns;
	m_funcRateSum.b = m_dLambdaIns+m_dLambdaDel*(m_pDeletionModel->MeanSize()-1.0);
	return true;
}

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