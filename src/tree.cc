// tree.cc - Copyright (C) 2004 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "tree.h"
#include "rand.h"

using namespace std;

////////////////////////////////////////////////////////////
//  class NewickNode
////////////////////////////////////////////////////////////

// Construct a NewickNode from the parser
NewickNode::NewickNode(NewickNode* p, const char *cs, double d) : m_dLen(d), m_pSub(p)
{
	if(cs)
		m_ssLabel = cs;
	else
		MakeName();
}

// Construct a name from the descendents of a node
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

Sequence::const_iterator Sequence::SeqPos(size_type uPos) const
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

Sequence::iterator Sequence::SeqPos(size_type uPos)
{
	iterator it = begin();
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

// Insert itBegin to itEnd at itPos

Sequence::size_type Sequence::Insertion(iterator itPos, const_iterator itBegin, const_iterator itEnd)
{
	if(itPos > end() || itPos < begin())
		return 0;

	size_type uRet = (size_type)(itEnd-itBegin);
	insert(itPos, itBegin, itEnd);
	m_uLength += uRet;
	return uRet;
}

// Delete uSize nucleotides at itBegin
Sequence::size_type Sequence::Deletion(iterator itBegin, size_type uSize)
{
	size_type uRet = 0;
	for(;uRet < uSize && itBegin != end(); ++itBegin)
	{
		// Skip Gaps
		if(itBegin->IsDeletion())
			continue;
		// Mark as Deleted-Root or Deleted-Insertion
		itBegin->SetType(itBegin->GetType()|Nucleotide::TypeDel);
		 ++uRet;
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
	ss.clear();
	for(const_iterator cit = begin(); cit != end(); ++cit)
		ss.push_back(cit->ToChar());
}

bool Nucleotide::FromChar(char ch)
{
	switch(ch&0xDF)
	{
	case 'A':
		m_ucNuc = NumAdenine;
		return true;
	case 'C':
		m_ucNuc = NumCytosine;
		return true;
	case 'G':
		m_ucNuc = NumGuanine;
		return true;
	case 'T':
		m_ucNuc = NumThymine;
		return true;
	}
	return false;
}
char Nucleotide::ToChar() const
{
	static const char csNuc[]	= "ACGT";
	static const char csType[]	= " +-=";
	return IsType(TypeRoot) ? csNuc[GetBase()] : csType[GetType() >> 2];
}


////////////////////////////////////////////////////////////
//  class Tree::Node
////////////////////////////////////////////////////////////

// Get the total sequence length of the node
Sequence::size_type Tree::Node::SeqLength() const
{
	Sequence::size_type uRet = 0;
	for(vector<Sequence>::const_iterator it = m_vSections.begin(); it != m_vSections.end(); ++it)
		uRet += it->SeqLength();
	return uRet;
}

// Flatten sections into one sequence
void Tree::Node::Flatten(Sequence& seq) const
{
	for(vector<Sequence>::const_iterator cit = m_vSections.begin();
		cit != m_vSections.end(); ++cit)
		seq.Append(*cit);
}

Tree::Node::iterator Tree::Node::SeqPos(Sequence::size_type uPos)
{
	vector<Sequence>::iterator itA;
	// Find section containing uPos
	for(itA = m_vSections.begin(); itA != m_vSections.end(); ++itA)
	{
		if(uPos < itA->SeqLength())
			break;
		uPos -= itA->SeqLength();
	}
	// Find actual iterator of uPos
	Sequence::iterator itB;
	if(itA != m_vSections.end())
		itB = itA->SeqPos(uPos);
	return iterator(itA,itB);
}

Tree::Node::const_iterator Tree::Node::SeqPos(Sequence::size_type uPos) const
{
	vector<Sequence>::const_iterator itA;
	// Find section containing uPos
	for(itA = m_vSections.begin(); itA != m_vSections.end(); ++itA)
	{
		if(uPos < itA->SeqLength())
			break;
		uPos -= itA->SeqLength();
	}
	// Find actual iterator of uPos
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
	// Construct the ur-root node if it doesn't exist
	m_map["_R()()T"];
	// process the newick tree beginning at its root
	ProcessNewickNode(pNode, "_R()()T");
	// increase number of sections
	m_nSec++;
}

void Tree::ProcessNewickNode(NewickNode* pNode, const string &ssAnc)
{
	// Process all sibs of the Newick Node
	if(pNode->m_pSib.get())
		ProcessNewickNode(pNode->m_pSib.get(), ssAnc);
	
	// Get a reference to this node and set up parameters
	Node& node = m_map[pNode->m_ssLabel];
	node.m_mBranchLens[ssAnc] = pNode->m_dLen;
	node.m_vAncestors.resize(m_nSec+1, ssAnc);
	//node.m_vAncestors[m_nSec] = ssAnc;
	
	// add node to tips if it is a tip and hasn't been accessed yet
	if(node.m_ssName.empty() && !pNode->m_pSub.get())
		m_vTips.push_back(pNode->m_ssLabel);
	
	// Give the node its name if it doesn't have one yet
	if(node.m_ssName.empty())
		node.m_ssName = pNode->m_ssLabel;
	
	// Process children
	if(pNode->m_pSub.get())
		ProcessNewickNode(pNode->m_pSub.get(), node.m_ssName);
}

// Evolve the sequences in the tree
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
	rNode.m_ssName = "_R()()T";
	rNode.m_bTouched = true;
	// Load sectiosn from template
	rNode.m_vSections = m_vDNASeqs;
	// Process Template
	for(vector<Sequence>::iterator it = rNode.m_vSections.begin();
		it != rNode.m_vSections.end(); ++it)
	{
		for(unsigned int u = 0;u<it->size();++u)
		{
			if((*it)[u].m_dRate < 0.0)
				(*it)[u].m_dRate = RandomRate(u);
			if((*it)[u].m_ucNuc >= 4)
				(*it)[u].m_ucNuc = RandomBase();
		}
	}
	// Evolve each tip
	for(vector<string>::const_iterator cit = m_vTips.begin(); cit != m_vTips.end(); ++cit)
		Evolve(m_map[*cit]);
}

// Evolve node rNode
void Tree::Evolve(Node &rNode)
{
	// do nothing if it has already been touched
	if(rNode.m_bTouched)
		return;
	rNode.m_bTouched = true;
	// Temporary Sequences
	map<string, Node> mapSeqs;
	// Evolve ancestors and assemble
	for(vector<string>::size_type a = 0; a < rNode.m_vAncestors.size(); ++a)
	{
		string &ssA = rNode.m_vAncestors[a];
		if(mapSeqs.find(ssA) == mapSeqs.end())
		{
			// Touch Ancestor, make sure it exists
			Node &aNode = m_map[ssA];
			Evolve(aNode);
			
			// Copy ancestor to temporary location
			mapSeqs[ssA] = aNode;
			// Evolve temporary location
			Evolve(mapSeqs[ssA], m_dTreeScale*rNode.m_mBranchLens[ssA]);
		}
		// Assemble final sequence
		rNode.m_vSections.push_back(mapSeqs[ssA].m_vSections[a]);
	}
}

// Evolve rNode a specific time 
void Tree::Evolve(Node &rNode, double dTime)
{
	dTime = fabs(dTime);
	if(dTime < DBL_EPSILON)
		return; // Nothing to evolve
	
	// Substitutions
	unsigned int uNuc = 0;
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
			// if dTemp is different from the previous one, recalculate probability matrix
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
			// get the base of the current nucleotide and pick new base
			unsigned int uBase = jt->GetBase();
			dTemp = rand_real();
			if(dTemp <= m_matSubst(uBase, 0))
				jt->SetBase(0);
			else if(dTemp <= m_matSubst(uBase, 1))
				jt->SetBase(1);
			else if(dTemp <= m_matSubst(uBase, 2))
				jt->SetBase(2);
			else
				jt->SetBase(3);

			++uNuc; // Increase position
		}
	}
	// Indel formation via Gillespie Algorithm

	// Check whether Indels are off
	if(m_dLambdaDel+m_dLambdaIns < DBL_EPSILON)
		return;

	// Get current length
	Sequence::size_type uLength = rNode.SeqLength()/m_uWidth;
	double dLength = (double)uLength;
	double dW = 1.0/m_funcRateSum(dLength);
	// Do indels
	for(double dt = rand_exp(dW); dt <= dTime; dt += rand_exp(dW))
	{
		// insertion or deletion
		if(rand_bool(m_funcRateIns(dLength)*dW))
		{
			//Insertion 
			Sequence::size_type ul = m_pInsertionModel->RandSize();
			Sequence::size_type uPos = (Sequence::size_type)rand_uint((uint32_t)uLength); // pos is in [0,L]
			// Construct sequence to be inserted
			Sequence seq;
			for(unsigned int uc = 0; uc < m_uWidth*ul; ++uc)
			{
				Nucleotide nuc = RandomNucleotide(uc);
				nuc.SetType(Nucleotide::TypeIns);
				seq.push_back(nuc);
			}
			// Find Position of Insertion
			Node::iterator itPos = rNode.SeqPos(uPos*m_uWidth);
			if(itPos.first == rNode.m_vSections.end())
			{
				// Insert at end of sequence
				uLength += rNode.m_vSections.back().Insertion(
					rNode.m_vSections.back().end(), seq.begin(), seq.end())/m_uWidth;
			}
			else
			{
				// Insert inside sequence
				uLength += itPos.first->Insertion(itPos.second, seq.begin(), seq.end())/m_uWidth;
			}
		}
		else if(uLength > 0)
		{
			// Deletion
			// Draw random size and random pos and rearrange
			Sequence::size_type ul = m_pDeletionModel->RandSize();
			Sequence::size_type uPos = rand_uint((uint32_t)(uLength+ul-2));
			Sequence::size_type uB = max(ul-1, uPos); 
			Sequence::size_type uSize = min(ul-1+uLength, uPos+ul)-uB;
			uB -= (ul-1);
			uB *= m_uWidth;
			uSize *= m_uWidth;
			// Find deletion point
			Node::iterator itPos = rNode.SeqPos(uB);
			Sequence::size_type uTemp = uSize;
			uTemp -= itPos.first->Deletion(itPos.second, uTemp);
			// Delete uSize nucleotides begin sensitive to gaps that overlap sections
			for(++itPos.first; uSize && itPos.first != rNode.m_vSections.end(); ++itPos.first)
				uTemp -= itPos.first->Deletion(itPos.first->begin(), uTemp);
			uLength -= (uSize-uTemp)/m_uWidth;
		}
		// update length
		dLength = (double)uLength;
		// new waiting time parameter
		dW = 1.0/m_funcRateSum(dLength);
	}
}

// Setup Evolutionary parameters
bool Tree::SetupEvolution(double pFreqs[], double pSubs[],
		const IndelModel::Params& rIns,
		const IndelModel::Params& rDel,
		unsigned int uWidth,
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
	
	// Setup Indel formation model
	// Insertion Rate
	m_dLambdaIns = rIns.dLambda;
	// Length Model
	if(m_dLambdaIns < DBL_EPSILON)
		m_pInsertionModel.release();
	else if(rIns.ssModel == "NB")
	{
		try {m_pInsertionModel.reset(new NegBnModel(rIns.vdModel));}
			catch(...) {return DawgError("Insertion model parameters not specified correctly.");}
	}
	else if(rIns.ssModel == "PL")
	{
		try {m_pInsertionModel.reset(new PowerModel(rIns.vdModel));}
			catch(...) {return DawgError("Insertion model parameters not specified correctly.");}
	}
	else if(rIns.ssModel == "US")
	{
		try {m_pInsertionModel.reset(new UserModel(rIns.vdModel));}
			catch(...) {return DawgError("Insertion model parameters not specified correctly.");}
	}
	else
		return DawgError("Unknown insertion model: \"%s\".", rDel.ssModel.c_str());

	// Deletion Rate
	m_dLambdaDel = rDel.dLambda;
	// Length model
	if(m_dLambdaDel < DBL_EPSILON)
		m_pDeletionModel.release();
	else if(rDel.ssModel == "NB")
	{
		try {m_pDeletionModel.reset(new NegBnModel(rDel.vdModel));}
			catch(...) {return DawgError("Deletion model parameters not specified correctly.");}
	}
	else if(rDel.ssModel == "PL")
	{
		try {m_pDeletionModel.reset(new PowerModel(rDel.vdModel));}
			catch(...) {return DawgError("Deletion model parameters not specified correctly.");}
	}
	else if(rDel.ssModel == "US")
	{
		try {m_pDeletionModel.reset(new UserModel(rDel.vdModel));}
			catch(...) {return DawgError("Deletion model parameters not specified correctly.");}
	}
	else
		return DawgError("Unknown deletion model: \"%s\".", rDel.ssModel.c_str());
	
	// Linear Functions for the Gillespie algorithm
	m_funcRateIns.m = m_dLambdaIns;
	m_funcRateIns.b = m_dLambdaIns;
	m_funcRateSum.m = m_dLambdaDel+m_dLambdaIns;
	m_funcRateSum.b = m_dLambdaIns+m_dLambdaDel*( (m_pDeletionModel.get()) ? m_pDeletionModel->MeanSize()-1.0 : 0.0);
	return true;
}

// Setup Root Template
bool Tree::SetupRoot(const std::vector<std::string> &vSeqs, const std::vector<unsigned int> &vLens,
					   const std::vector<std::vector<double> > &vRates)
{
	// Clear Template
	m_vDNASeqs.clear();
	
	// Check to see if sequence is specified
	if(vSeqs.size())
	{	
		// Read sequence of each section
		for(vector<string>::const_iterator cit = vSeqs.begin(); cit != vSeqs.end(); ++cit)
		{
			Sequence seq(BlockTrim((unsigned int)cit->size()));
			for(unsigned int u=0; u<seq.size(); ++u)
				if(!seq[u].FromChar((*cit)[u]))
					return DawgError("Unknown character, \"%c\", in Sequence", (*cit)[u]);
			m_vDNASeqs.push_back(seq);
		}
	}
	else
	{
		// Create random sequences
		for(vector<unsigned int>::const_iterator cit = vLens.begin(); cit != vLens.end(); ++cit)
			m_vDNASeqs.push_back(Sequence(*cit*m_uWidth));
	}
	// Check to see if rates are specified
	if(vRates.size())
	{
		// Read rates of each section
		for(unsigned int u=0; u < m_vDNASeqs.size(); ++u)
		{
			double dTemp = 0.0;
			for(unsigned int v=0; v < m_vDNASeqs[u].size(); ++v)
			{
				m_vDNASeqs[u][v].m_dRate = vRates[u][v];
				dTemp += vRates[u][v];

			}
			// Scale the Expected Rate to 1.0
			for(unsigned int v=0; v < m_vDNASeqs[u].size(); ++v)
				m_vDNASeqs[u][v].m_dRate /= dTemp;
		}
	}
	return true;
}

unsigned char Tree::RandomBase() const
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

double Tree::RandomRate(Sequence::size_type uPos) const
{
	uPos %= m_uWidth;
	if(m_vdIota[uPos] > DBL_EPSILON && rand_bool(m_vdIota[uPos]))
		return 0.0;  // Site Invariant
	else if(m_vdGamma[uPos] > DBL_EPSILON)
		return rand_gamma1(m_vdGamma[uPos]); // Gamma with mean 1.0 and var of m_dGamma
	else
		return 1.0;
}

void Tree::Align(Alignment &aln) const
{
	// construct a table of flattened sequences
	vector<Sequence> vTable;
	vector<string> vNames;
	for(Node::Map::const_iterator cit = m_map.begin(); cit != m_map.end(); ++cit)
	{
		vNames.push_back(cit->second.m_ssName);
		Sequence s;
		cit->second.Flatten(s);
		vTable.push_back(s);
	}
	// Alignment rules:
	// Insertion & Deleted Insertion  : w/ ins, deleted ins, or gap
	// Deletion & Original Nucleotide : w/ del, original nucl
	
	// States: Quit (0), Del/Root (1), Ins (2), InsDel (3)
	unsigned char uState = 1;
	// Go through each column, adding gaps where neccessary
	for(unsigned int uCol = 0; uState; uCol++)
	{
		// Set to quit
		uState = 0;
		for(vector<Sequence>::const_iterator cit = vTable.begin();
			cit != vTable.end(); ++cit)
		{
			if(uCol >= cit->size())
				continue;
			if(uState == 0)
				uState = 1;	// Nucleotide exists clear quit
			if((*cit)[uCol].IsType(Nucleotide::TypeIns))
			{
				// Gaps need to be added mark and break
				uState = 2;
				break;
			}
			else if((*cit)[uCol].IsType(Nucleotide::TypeDelIns))
			{
				uState = 3;
			}
		}
		if(uState == 3)
		{
			for(vector<Sequence>::iterator it = vTable.begin();
				it != vTable.end(); ++it)
			{
				if(uCol < it->size() && (*it)[uCol].GetType() == Nucleotide::TypeDelIns)
					it->erase(it->begin()+uCol);
			}
			uCol--;
		}
		else if(uState == 2)
		{
			// Add gaps where neccessary
			for(vector<Sequence>::iterator it = vTable.begin();
				it != vTable.end(); ++it)
			{
				if(uCol < it->size())
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
				else
					// Add gap to end of sequence
					it->resize(uCol+1, Nucleotide(Nucleotide::TypeIns, 1.0));
			}
		}
	}

	// Add aligned sequences to alingment set
	for(unsigned int u = 0; u < vNames.size(); ++u)
	{
		// Skip any sequence that begin with one of the two special characters
		if(vNames[u][0] == '(' || vNames[u][0] == '_')
			continue;
		// Add Sequence to alignment
		string& ss = aln[vNames[u]];
		vTable[u].ToString(ss);
	}
}


