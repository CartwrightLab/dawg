#include "subst.h"
#include "node.h"
#include "dawg.h"
#include "rand.h"

#include <math.h>
#include <float.h>

bool SubstProcessor::Setup(double pFreqs[], double pSubs[])
{
	if(pFreqs[0] < 0 || pFreqs[1] < 0 || pFreqs[2] < 0 || pFreqs[3] < 0)
		return DawgError("Nucleotide frequences need to be positive.");
	memcpy(m_dFreqs, pFreqs, 4*sizeof(double));
	m_dNucCumFreqs[0] = pFreqs[0];
	m_dNucCumFreqs[1] = pFreqs[1]+m_dNucCumFreqs[0];
	m_dNucCumFreqs[2] = pFreqs[2]+m_dNucCumFreqs[1];
	m_dNucCumFreqs[3] = 1.0;
	if(m_dNucCumFreqs[0] > 1.0 || m_dNucCumFreqs[1] > 1.0 || m_dNucCumFreqs[2] > 1.0)
		return DawgError("Nucleotide frequences need to sum to 1.");
	
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
	matQ.Scale(matQ, -1.0/(vecF[0]*matQ(0,0)+vecF[1]*matQ(1,1)+
		vecF[2]*matQ(2,2)+vecF[3]*matQ(3,3)));
	
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

void SubstProcessor::Process(Node *pNode)
{
	double dLen = pNode->ScaledLength();
	double dTemp;
	for(Seq::iterator nit = pNode->Sequence().begin();
		nit != pNode->Sequence().end(); ++nit)
	{
		if(nit->m_dRate < DBL_EPSILON)
			continue; // Invariant Site
		SetTime(dLen*nit->m_dRate);
		dTemp = rand_real();
		if(dTemp <= m_matSubst(nit->m_nuc, 0))
			nit->m_nuc = 0;
		else if(dTemp <= m_matSubst(nit->m_nuc, 1))
			nit->m_nuc = 1;
		else if(dTemp <= m_matSubst(nit->m_nuc, 2))
			nit->m_nuc = 2;
		else
			nit->m_nuc = 3;
	}
}

void SubstProcessor::SetTime(double dTime)
{
	if(dTime == m_dOldTime)
		return;
	m_dOldTime = dTime;
	Vector4  vec;
	vec[0] = exp(dTime*m_vecL[0]);
	vec[1] = exp(dTime*m_vecL[1]);
	vec[2] = exp(dTime*m_vecL[2]);
	vec[3] = exp(dTime*m_vecL[3]);
	Matrix44 mat; mat.Scale(vec, m_matU);
	m_matSubst.Multiply(m_matV, mat);
	for(Matrix44::Pos i=0;i<4;++i)
	{
		m_matSubst(i,1) += m_matSubst(i,0);
		m_matSubst(i,2) += m_matSubst(i,1);
		//m_matSubst(i,3) = 1.0;
	}
}