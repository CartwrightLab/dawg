#include "indel.h"
#include "dawg.h"
#include "rand.h"

#include <float.h>
#include <math.h>

using namespace std;

// rate of del: lambda*length+lambda*(E(size)-1)
// rate of ins: lambda*length+lambda

////////////////////////////////////////////////////////////
//  class IndelProcessor
////////////////////////////////////////////////////////////

bool IndelProcessor::Setup(const IndelModel::Params& rIns, const IndelModel::Params& rDel)
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

void IndelProcessor::Process(Node *pNode)
{
	// Check to see if we are doing indels
	if(m_funcRateSum(0.0) < DBL_EPSILON)
		return;

	double dtMax = fabs(pNode->ScaledLength());
	unsigned long uSize = pNode->Sequence().size();
	double dSize = (double)uSize;
	double dW = 1.0/m_funcRateSum(dSize);
	double dt = rand_exp(dW);
	while(dt <= dtMax)
	{
		if(rand_bool(m_funcRateIns(dSize)*dW))
		{
			//Insertion
			pNode->Insert(rand_ulong(uSize), m_pInsertionModel->RandSize());
		}
		else
		{
			//Deletion
			unsigned long ul = m_pDeletionModel->RandSize();
			pNode->Delete(rand_ulong(uSize+ul-1), ul);
		}
		uSize = pNode->Sequence().size();
		dSize = (double)uSize;
		dW = 1.0/m_funcRateSum(dSize);
		dt += rand_exp(dW);
	}
}

////////////////////////////////////////////////////////////
//  class NegBnModel
////////////////////////////////////////////////////////////

NegBnModel::NegBnModel(const vector<double>& vdModel)
{
	if(vdModel.size() < 2)
		throw(DawgError("Negative Binomial Model requires two parameters."));
	m_uR = (unsigned int)vdModel[0];
	m_dQ = vdModel[1];
	if(m_uR == 0)
		throw(DawgError("Negative Binomial Model requires R > 0."));
	if(0.0 > m_dQ || m_dQ > 1.0)
		throw(DawgError("Negative Binomial Model requires that Q be a a probability."));
}
unsigned long NegBnModel::RandSize() const
{
	return 1ul+rand_negbinomial(m_uR, m_dQ);
}
double NegBnModel::MeanSize() const
{
	return 1.0+m_uR*m_dQ/(1.0-m_dQ);
}

////////////////////////////////////////////////////////////
//  class UserModel
////////////////////////////////////////////////////////////

UserModel::UserModel(const vector<double>& vdModel)
{
	m_vSizesCum.clear();
	double dSum = 0.0;
	for(vector<double>::const_iterator cit = vdModel.begin();
		cit != vdModel.end(); ++cit)
	{
		if(*cit > 1.0 || dSum > 1.0)
			throw(DawgError("User Gap Model parameters must sum to one."));
		dSum += *cit;
		m_vSizesCum.push_back(dSum);
	}
}
unsigned long UserModel::RandSize() const
{
	double d = rand_real();
	unsigned long ul = 1ul;
	for(vector<double>::const_iterator dit = m_vSizesCum.begin();
		dit != m_vSizesCum.end() && d > *dit; ++dit)
		++ul;
	return ul;
}
double UserModel::MeanSize() const
{
	double dSize = 1.0;
	double dLast = 0.0;
	double dSum = 0.0;
	for(vector<double>::const_iterator dit = m_vSizesCum.begin();
		dit != m_vSizesCum.end(); ++dit)
	{
		dSum += dSize*(*dit-dLast);
		dLast = *dit;
		dSize += 1.0;
	}
	return dSum;
}