// indel.cc - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

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
IndelModel::size_type NegBnModel::RandSize() const
{
	return (IndelModel::size_type)(1ul+rand_negbinomial(m_uR, m_dQ));
}
double NegBnModel::MeanSize() const
{
	return 1.0+m_uR*m_dQ/(1.0-m_dQ);
}

////////////////////////////////////////////////////////////
//  class UserModel
////////////////////////////////////////////////////////////

UserModel::UserModel(const vector<double>& vdModel) : m_dMean(0.0)
{
	m_vSizesCum.clear();
	double dSum = 0.0, dSize = 1.0;
	
	for(vector<double>::const_iterator cit = vdModel.begin();
		cit != vdModel.end(); ++cit)
	{
		if(*cit < 0.0 || *cit > 1.0 || dSum > 1.0)
			throw(DawgError("User Model parameters must be positive and sum to one."));
		m_dMean += *cit*dSize;
		dSum += *cit;
		m_vSizesCum.push_back(dSum);
		dSize += 1.0;
	}
}

IndelModel::size_type UserModel::RandSize() const
{
	double d = rand_real();
	IndelModel::size_type ul = 1ul;
	
	// not efficient, but it gets the job done.
	for(vector<double>::const_iterator dit = m_vSizesCum.begin();
		dit != m_vSizesCum.end() && d >= *dit; ++dit)
		++ul;
	
	return ul;
}
double UserModel::MeanSize() const
{
	return m_dMean;
}

////////////////////////////////////////////////////////////
//  class PowerModel
////////////////////////////////////////////////////////////
PowerModel::PowerModel(const vector<double>& vdModel) : m_dMean(0.0)
{
	if(vdModel.size() < 2)
		throw(DawgError("Power Law Model requires two parameters."));

	m_dA = vdModel[0];
	m_uM = (IndelModel::size_type)vdModel[1];
	
	if(m_dA < 1.0)
		throw(DawgError("Power Law Model requires a >= 1.0."));
	if(vdModel[1] < 1.0)
		throw(DawgError("Power Law Model requires Max >= 1."));
	
	double dSum = 0.0, d;
	for(IndelModel::size_type u = 1;u<=m_uM;++u)
	{
		d = pow((double)u, -m_dA);
		m_dMean += u*d;
		dSum += d;
	}
	m_dMean /= dSum;
}

IndelModel::size_type PowerModel::RandSize() const
{
	IndelModel::size_type u;
	do {
		u = (IndelModel::size_type)rand_zipf(m_dA);
	} while (u > m_uM);
	return u;
}

double PowerModel::MeanSize() const
{
	return m_dMean;
}

