#include "dawg.h"
#include "indel.h"
#include "rand.h"

using namespace std;

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