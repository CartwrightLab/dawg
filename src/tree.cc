// tree.cc - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

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

Sequence::const_iterator Sequence::SeqPos(unsigned long uPos) const
{
	const_iterator it = begin();
	// Skip deletions
	while(it->IsDeletion()) {++it;}
	while(uPos--)
	{
		++it;
		// Skip deletions
		while(it->IsDeletion()) {++it;}
	}
	return it;
}

Sequence::iterator Sequence::SeqPos(unsigned long uPos)
{
	iterator it = begin();
	// Skip deletions in the history
	while(it->IsDeletion()) {++it;}
	while(uPos--)
	{
		++it;
		// Skip deletions in the history
		while(it->IsDeletion()) {++it;}
	}
	return it;
}

unsigned long Sequence::Insertion(iterator itPos, const_iterator itBegin, const_iterator itEnd)
{
	unsigned long uRet = (unsigned long)(itEnd-itBegin);
	size_type off = size() == 0 ? 0 : itPos - begin();
	insert(itPos, itBegin, itEnd);
	iterator it = begin()+off;
	iterator itE = it+uRet;
	for(;it!=itE;++it)
		it->SetType(Nucleotide::TypeIns);
	m_uLength += uRet;
	return uRet;
}

unsigned long Sequence::Deletion(iterator itBegin, unsigned long uSize)
{
	unsigned long uRet = 0;
	for(;uRet < uSize && itBegin != end(); ++itBegin)
	{
		if(itBegin->IsDeletion())
			continue;
		itBegin->SetType(itBegin->IsInsertion() ?
			Nucleotide::TypeDelIns : Nucleotide::TypeDel);
		uRet++;
	}
	m_uLength -= uRet;
	return uRet;
}


void Sequence::Append(const Sequence &seq)
{
		insert(end(), seq.begin(), seq.end());
		m_uLength += seq.m_uLength;
}

void Sequence::ToString(std::string &ss) const
{
	for(const_iterator cit = begin(); cit != end(); ++cit)
		ss.push_back(cit->ToChar());
}


////////////////////////////////////////////////////////////
//  class Tree::Node
////////////////////////////////////////////////////////////

unsigned long Tree::Node::SeqLength() const
{
	unsigned long uRet = 0;
	for(vector<Sequence>::const_iterator it = m_vSections.begin(); it != m_vSections.end(); ++it)
		uRet += it->SeqLength();
	return uRet;
}

void Tree::Node::Flatten(Sequence& seq) const
{
	for(vector<Sequence>::const_iterator cit = m_vSections.begin();
		cit != m_vSections.end(); ++cit)
		seq.Append(*cit);
}

Tree::Node::iterator Tree::Node::SeqPos(unsigned long uPos)
{
	vector<Sequence>::iterator itA;
	for(itA = m_vSections.begin(); itA != m_vSections.end(); ++itA)
	{
		if(uPos < itA->SeqLength())
			break;
		uPos -= itA->SeqLength();
	}
	Sequence::iterator itB;
	if(itA != m_vSections.end())
		itB = itA->SeqPos(uPos);
	return iterator(itA,itB);
}

Tree::Node::const_iterator Tree::Node::SeqPos(unsigned long uPos) const
{
	vector<Sequence>::const_iterator itA;
	for(itA = m_vSections.begin(); itA != m_vSections.end(); ++itA)
	{
		if(uPos < itA->SeqLength())
			break;
		uPos -= itA->SeqLength();
	}
	Sequence::const_iterator itB;
	if(itA != m_vSections.end())
		itB = itA->SeqPos(uPos);
	return const_iterator(itA,itB);
}


////////////////////////////////////////////////////////////
//  class Tree
////////////////////////////////////////////////////////////

void Tree::ProcessTree(NewickNode* pNode)
{
	m_map["_R()()T"];
	ProcessNewickNode(pNode, m_map.find("_R()()T"));
	m_nSec++;
}

void Tree::ProcessNewickNode(NewickNode* pNode, Node::Handle hAnc)
{
	if(pNode->m_pSib.get())
		ProcessNewickNode(pNode->m_pSib.get(), hAnc);

	Node& node = m_map[pNode->m_ssLabel];
	node.m_mBranchLens[hAnc] = pNode->m_dLen;
	node.m_vAncestors.resize(m_nSec+1, m_map.end());
	node.m_vAncestors[m_nSec] = hAnc;
			
	if(pNode->m_pSub.get())
		ProcessNewickNode(pNode->m_pSub.get(), m_map.find(pNode->m_ssLabel));
}

void Tree::Evolve()
{
	// Reset Sequences
	for(Node::Map::iterator it=m_map.begin(); it!=m_map.end();++it)
	{
		it->second.m_vSections.clear();
		it->second.m_bTouched = false;
	}
	// Setup Root
	Node& rNode = m_map["_R()()T"];
	rNode.m_bTouched = true;
	rNode.m_vSections = m_vDNASeqs;
	for(vector<Sequence>::iterator it = rNode.m_vSections.begin();
		it != rNode.m_vSections.end(); ++it)
	{
		for(unsigned int u = 0;u<it->size();++u)
		{
			if((*it)[u].m_dRate < 0.0)
				(*it)[u].m_dRate = RandomRate(u);
			if((*it)[u].m_ucNuc >= 4)
				(*it)[u].m_ucNuc = RandomNuc();
		}
	}
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
		if(mapSeqs.find(rNode.m_vAncestors[a]) == mapSeqs.end())
		{
			// Touch Ancestor, make sure it exists
			Evolve(rNode.m_vAncestors[a]->second);
			// Copy ancestor to temporary location
			mapSeqs[rNode.m_vAncestors[a]] = rNode.m_vAncestors[a]->second;
			// Evolve temporary location
			Evolve(mapSeqs[rNode.m_vAncestors[a]], m_dTreeScale*rNode.m_mBranchLens[rNode.m_vAncestors[a]]);
		}
		// Assemble final sequence
		rNode.m_vSections.push_back(mapSeqs[rNode.m_vAncestors[a]].m_vSections[a]);
	}
}

void Tree::Evolve(Node &rNode, double dTime)
{
	dTime = fabs(dTime);
	if(dTime < DBL_EPSILON)
		return; // Nothing to evolve
	
	// Substitutions
	unsigned long uNuc = 0;
	for(vector<Sequence>::iterator it = rNode.m_vSections.begin(); it != rNode.m_vSections.end(); ++it)
	{
		for(Sequence::iterator jt = it->begin(); jt != it->end(); ++jt)
		{
			// Skip any position that is a deletion
			if(jt->IsDeletion())
				continue;
			// Total Evolution Rate for the position
			double dTemp = dTime*jt->m_dRate*m_vdScale[uNuc%m_uWidth];
			if(dTemp < DBL_EPSILON)
				continue; // Invariant Site
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
			if(dTemp <= m_matSubst(jt->m_ucNuc, 0))
				jt->m_ucNuc = 0;
			else if(dTemp <= m_matSubst(jt->m_ucNuc, 1))
				jt->m_ucNuc = 1;
			else if(dTemp <= m_matSubst(jt->m_ucNuc, 2))
				jt->m_ucNuc = 2;
			else
				jt->m_ucNuc = 3;

			++uNuc; // Increase position
		}
	}
	// Indels
	if(m_funcRateSum(0.0) < DBL_EPSILON)
		return;

	unsigned long uLength = rNode.SeqLength()/m_uWidth;
	double dLength = (double)uLength;
	double dW = 1.0/m_funcRateSum(dLength);
	for(double dt = rand_exp(dW); dt <= dTime; dt += rand_exp(dW))
	{
		if(rand_bool(m_funcRateIns(dLength)*dW))
		{
			//Insertion
			unsigned long ul = m_pInsertionModel->RandSize();
			unsigned long uPos = rand_ulong(uLength); // pos is in [0,L]
			Sequence seq;
			for(unsigned int uc = 0; uc < m_uWidth*ul; ++uc)
				seq.push_back(RandomNucleotide(uc));
			Node::iterator itPos = rNode.SeqPos(uPos*m_uWidth);
			if(itPos.first == rNode.m_vSections.end())
			{
				uLength += rNode.m_vSections.back().Insertion(
					rNode.m_vSections.back().end(), seq.begin(), seq.end())/m_uWidth;
			}
			else
			{
				uLength += itPos.first->Insertion(itPos.second, seq.begin(), seq.end())/m_uWidth;
			}
		}
		else
		{
			//Deletion
			unsigned long ul = m_pDeletionModel->RandSize();
			unsigned long uPos = rand_ulong(uLength+ul-1)+1;
			unsigned long uB = (ul >= uPos) ? 0 : uPos - ul;
			unsigned long uSize = (uPos > uLength) ? uLength : uPos;
			uSize -= uB;
			uSize *= m_uWidth;
			uB *= m_uWidth;

			Node::iterator itPos = rNode.SeqPos(uB);
			uSize -= itPos.first->Deletion(itPos.second, uSize);
			for(++itPos.first; uSize && itPos.first != rNode.m_vSections.end(); ++itPos.first)
				uSize -= itPos.first->Deletion(itPos.first->begin(), uSize);
			uLength -= (ul*m_uWidth - uSize)/m_uWidth;
		}
		dLength = (double)uLength;
		dW = 1.0/m_funcRateSum(dLength);
	}
}

bool Tree::SetupEvolution(double pFreqs[], double pSubs[],
		const IndelModel::Params& rIns,
		const IndelModel::Params& rDel,
		unsigned long uWidth,
		const std::vector<double> &vdGamma,
		const std::vector<double> &vdIota,
		const std::vector<double> &vdScale,
		double dTreeScale)
{
	// Verifiy Parameters
	if(pFreqs[0] < 0.0 || pFreqs[1] < 0.0 || pFreqs[2] < 0.0 || pFreqs[3] < 0.0)
		return DawgError("Nucleotide frequences need to be positive.");
	pFreqs[3] = 1.0-pFreqs[0]-pFreqs[1]-pFreqs[2];
	if( pFreqs[3] < 0.0 )
		return DawgError("Nucleotide frequencies need to sum to 1.0.");
	if(pSubs[0] < 0.0 || pSubs[1] < 0.0 || pSubs[2] < 0.0
		|| pSubs[3] < 0.0 || pSubs[4] < 0.0 || pSubs[5] < 0.0)
		return DawgError("Substitution rates need to be positive.");

	if(rIns.dLambda < 0.0)
		return DawgError("Lambda (Ins) must not be negative.");
	if(rDel.dLambda < 0.0)
		return DawgError("Lambda (Del) must not be negative.");
	if(uWidth == 0)
		return DawgError("Width must be positive.");
	if(vdGamma.size() != uWidth)
		return DawgError("Gamma must have the same size as the value of Width.");
	if(vdIota.size() != uWidth)
		return DawgError("Iota must have the same size as the value of Width.");
	if(vdScale.size() != uWidth)
		return DawgError("Scale must have the same size as the value of Width.");
	for(vector<double>::const_iterator cit = vdGamma.begin(); cit != vdGamma.end(); ++cit)
	{
		if(*cit < 0.0)
			return DawgError("Invalid Gamma, \"%f\".  Gamma must be positive.", *cit);
	}
	for(vector<double>::const_iterator cit = vdIota.begin(); cit != vdIota.end(); ++cit)
	{
		if(0.0 > *cit || *cit > 1.0)
			return DawgError("Invalid Iota, \"%f\".  Iota must be a probability.", *cit);
	}
	for(vector<double>::const_iterator cit = vdScale.begin(); cit != vdScale.end(); ++cit)
	{
		if(*cit <= 0.0)
			return DawgError("Invalid Scale, \"%f\". Scale must be positive.", *cit);
	}
	if(dTreeScale <= 0.0)
		return DawgError("Invalid TreeScale, \"%f\". TreeScale must be positive.", dTreeScale);

	// Setup Frame
	m_uWidth = uWidth;

	// Setup Rate Parameters
	m_vdGamma = vdGamma;
	m_vdIota = vdIota;

	// Setup Scale
	m_vdScale = vdScale;

	// Setup TreeScale
	m_dTreeScale = dTreeScale;

	// Setup Cumulative Frequencies
	m_dNucCumFreqs[0] = pFreqs[0];
	m_dNucCumFreqs[1] = m_dNucCumFreqs[0]+pFreqs[1];
	m_dNucCumFreqs[2] = m_dNucCumFreqs[1]+pFreqs[2];
	m_dNucCumFreqs[3] = 1.0;

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
	double dX = 0.0;
	for(unsigned int i=0;i<m_vdIota.size();++i)
		dX -= (1.0-m_vdIota[i])*m_vdScale[i];
	dX = m_vdIota.size()/dX;
	matQ(0,0) = -(matQ(0,1)+matQ(0,2)+matQ(0,3));
	matQ(1,1) = -(matQ(1,0)+matQ(1,2)+matQ(1,3));
	matQ(2,2) = -(matQ(2,0)+matQ(2,1)+matQ(2,3));
	matQ(3,3) = -(matQ(3,0)+matQ(3,1)+matQ(3,2));
	matQ.Scale(matQ, dX/((vecF[0]*matQ(0,0)+vecF[1]*matQ(1,1)+
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

	m_dLambdaIns = rIns.dLambda;
	if(m_dLambdaIns < DBL_EPSILON)
		m_pInsertionModel.release();
	else if(rIns.ssModel == "NB")
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
	if(m_dLambdaDel < DBL_EPSILON)
		m_pDeletionModel.release();
	else if(rDel.ssModel == "NB")
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
	m_funcRateSum.b = m_dLambdaIns+m_dLambdaDel*( (m_pDeletionModel.get()) ? m_pDeletionModel->MeanSize()-1.0 : 0.0);
	return true;
}

bool Tree::SetupRoot(const std::vector<std::string> &vSeqs, const std::vector<unsigned long> &vLens,
					   const std::vector<std::vector<double> > &vRates)
{
	m_vDNASeqs.clear();
	if(vSeqs.size())
	{	
		for(vector<string>::const_iterator cit = vSeqs.begin(); cit != vSeqs.end(); ++cit)
		{
			Sequence seq(BlockTrim((unsigned long)cit->size()));
			for(unsigned int u=0; u<seq.size(); ++u)
				if(!seq[u].FromChar((*cit)[u]))
					return DawgError("Unknown character, \"%c\", in Sequence", (*cit)[u]);
			m_vDNASeqs.push_back(seq);
		}
	}
	else
	{
		for(vector<unsigned long>::const_iterator cit = vLens.begin(); cit != vLens.end(); ++cit)
			m_vDNASeqs.push_back(Sequence(BlockTrim(*cit)));
	}
	if(vRates.size())
	{
		for(unsigned int u=0; u < m_vDNASeqs.size(); ++u)
			for(unsigned int v=0; v < m_vDNASeqs[u].size(); ++v)
				m_vDNASeqs[u][v].m_dRate = vRates[u][v];
	}
	return true;
}

unsigned char Tree::RandomNuc() const
{
	double d = rand_real();
	if(d <= m_dNucCumFreqs[0])
		return 0; // A
	else if(d <= m_dNucCumFreqs[1])
		return 1; // C
	else if(d <= m_dNucCumFreqs[2])
		return 2; // G
	else
		return 3; // T
}

double Tree::RandomRate(unsigned long uPos) const
{
	uPos %= m_uWidth;
	if(m_vdIota[uPos] > DBL_EPSILON && rand_bool(m_vdIota[uPos]))
		return 0.0;  // Site Invariant
	else if(m_vdGamma[uPos] > DBL_EPSILON)
		return rand_gamma1(m_vdGamma[uPos]); // Gamma with mean 1.0 and var of m_dGamma
	else
		return 1.0;
}

void Tree::Align(Alignment &aln, bool bGapPlus, bool bGapSingleChar, bool bLowerCase) const
{
	// construct a table of flattened sequences
	vector<Sequence> vTable;
	vector<string> vNames;
	for(Node::Map::const_iterator cit = m_map.begin(); cit != m_map.end(); ++cit)
	{
		vNames.push_back(cit->first);
		Sequence s;
		cit->second.Flatten(s);
		vTable.push_back(s);
	}
	// Alignment rules:
	// Insertion & Deleted Insertion  : w/ ins, deleted ins, or gap
	// Deletion & Original Nucleotide : w/ del, original nucl
	
	unsigned char uState = 1;
	for(unsigned int uCol = 0; uState; uCol++)
	{
		uState = 0;
		for(vector<Sequence>::const_iterator cit = vTable.begin();
			cit != vTable.end(); ++cit)
		{
			if(uCol >= cit->size())
				continue;
			uState = 1;
			if((*cit)[uCol].IsInsertion())
			{
				uState = 2;
				break;
			}
		}
		if(uState == 2)
		{
			for(vector<Sequence>::iterator it = vTable.begin();
				it != vTable.end(); ++it)
			{
				switch((*it)[uCol].GetType())
				{
				case Nucleotide::TypeIns:
					(*it)[uCol].SetType(Nucleotide::TypeRoot);
					break;
				case Nucleotide::TypeDelIns:
					break;
				default:
					it->insert(it->begin()+uCol, Nucleotide(Nucleotide::TypeIns, 1.0));
					break;
				}
			}
		}


	}
	for(unsigned int u = 0; u < vNames.size(); ++u)
	{
		string& ss = aln[vNames[u]];
		vTable[u].ToString(ss);
		//GapPlus and GappSingleChar
		int nGapState = 0;
		for(unsigned int v = 0; v < ss.length(); ++v)
		{
			if(bGapSingleChar)
			{
				if(ss[v] == '-' || ss[v] == '=')
				{	
					if(nGapState == 1)
						ss[v] = '?';
					else
						nGapState = 1;
				}
				else if(ss[v] == '+')
				{
					if(nGapState == 2)
						ss[v] = '?';
					else
						nGapState = 2;
				}
				else
					nGapState = 0;
			}
			if(!bGapPlus && (ss[v] == '+' || ss[v] == '='))
				ss[v] = '-';
			if(bLowerCase && 0x41 <= ss[v] && ss[v] <= 0x5A)
				ss[v] |= 0x20;
		}
	}

}


