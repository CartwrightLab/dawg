#pragma warning(disable: 4512)

#include <algorithm>
#include "indelmodel.h"

bool Dawg::IndelDistribution::Create(Dawg::IndelDistribution::Type mod, const std::vector<double> &vdParams)
{
	m_vdParams = vdParams;
	m_model = mod;
	switch(m_model)
	{
	case MOD_USER:
		rand_user.distribution() = discrete_distribution<result_type>(vdParams.begin(), vdParams.end());
		return true;
	case MOD_GEO:
		rand_geo.distribution() = boost::geometric_distribution<result_type>(vdParams[0]);
		return true;
	case MOD_ZIPF:
		rand_zipf.distribution() = truncated_zipf_distribution<result_type>(vdParams[0], result_type(vdParams[1]));
		return true;
	default:
		return false;
	};
}

Dawg::IndelDistribution::result_type Dawg::IndelDistribution::Rand()
{
	switch(m_model)
	{
	case MOD_USER:
		return rand_user();
	case MOD_GEO:
		return rand_geo();
	case MOD_ZIPF:
	default:
		return rand_zipf();
	};
}

bool Dawg::IndelModel::Create(double dInsRate, IndelDistribution::Type modIns, const std::vector<double> &vdInsParams,
				double dDelRate, IndelDistribution::Type modDel, const std::vector<double> &vdDelParams,
				const SequenceFactory &seqfac)
{
	m_dInsRate = dInsRate;
	m_dDelRate = dDelRate;
	if(!rand_ins.Create(modIns, vdInsParams) ||
	   !rand_del.Create(modDel, vdDelParams))
	   return false;
	//rand_seq = seqfac;
	return true;
}

void Dawg::IndelModel::Process(Dawg::Sequence &seq, double dTime)
{
	if(m_dInsRate+m_dDelRate < DBL_EPSILON+DBL_EPSILON || dTime < DBL_EPSILON)
		return;
	double dLength = (double)seq.Length();
	const double dM = m_dInsRate+m_dDelRate;
	const double dB = m_dInsRate;
	double dW = dLength*dM+dB;
	double dP = (m_dInsRate*dLength+m_dInsRate)/dW;
	Sequence seqIns;
	for(double dT = rand_exp(dW); dT <= dTime; dT += rand_exp(dW))
	{
		if(rand_bool(dP))
		{
			//do Insertion
			rand_seq(seqIns, rand_ins());
			seq.Insert(rand_uint(seq.Length()), seqIns);
		}
		else
		{
			//do Deletion
			Sequence::pos_type uLen, uP, uB, uE;
			uLen = rand_del();
			if(uLen == Sequence::pos_type(1))
			{
				// accelerate most common type
				uB = rand_uint(seq.Length()-1);
				uE = uB+1;
			}
			else
			{
				uP = rand_uint(seq.Length()+uLen-2)+1;
				uB = std::max(uLen, uP)-uLen; 
				uE = std::min(seq.Length(), uP);
			}
			seq.Delete(uB, uE);
		}
		
		//Update Gillisepe Algorithm
		dLength = (double)seq.Length();
		dW = dLength*dM+dB;
		dP = (m_dInsRate*dLength+m_dInsRate)/dW;
	}
}

