#ifndef DAWG_SUBST_H
#define DAWG_SUBST_H

#include "matrix.h"

class Node;

//Substitution Model is GTR+I+G
class SubstProcessor
{
public:
	SubstProcessor() : m_dOldTime(-1.0) { }
	bool Setup(double pFreqs[], double pSubs[]);
	void Process(Node *pNode);
	void SetTime(double dTime);

protected:
	double m_dOldTime;
	double m_dFreqs[4];
	double m_dNucCumFreqs[4];
	Matrix44 m_matSubst;
	Matrix44 m_matV;
	Matrix44 m_matU;
	Matrix44 m_matQ;
	Matrix44 m_matR;
	Vector4  m_vecL;
};

#endif