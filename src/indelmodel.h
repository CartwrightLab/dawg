#ifndef DAWG_INDELMODEL_H
#define DAWG_INDELMODEL_H

#include "rand.h"
#include "sequence.h"

namespace Dawg {

class IndelDistribution
{
public:
	typedef unsigned int result_type;
	enum Type {MOD_ZIPF, MOD_GEO, MOD_USER};
	
	IndelDistribution() : m_model(MOD_ZIPF),
		rand_zipf(g_rng, truncated_zipf_distribution<result_type>()),
		rand_user(g_rng, discrete_distribution<result_type>()),
		rand_geo(g_rng, boost::geometric_distribution<result_type>())
	{ }

	bool Create(Type mod, const std::vector<double> &vdParams);

	result_type Rand();
	result_type operator()() { return Rand(); }

private:
	Type m_model;
	std::vector<double> m_vdParams;

	boost::variate_generator<DawgRng&, truncated_zipf_distribution<result_type> > rand_zipf;
	boost::variate_generator<DawgRng&, discrete_distribution<result_type> > rand_user;
	boost::variate_generator<DawgRng&, boost::geometric_distribution<result_type> > rand_geo;

};

class IndelModel
{
public:
	IndelModel() { }

	bool Create(double dInsRate, IndelDistribution::Type modIns, const std::vector<double> &vdInsParams,
				double dDelRate, IndelDistribution::Type modDel, const std::vector<double> &vdDelParams,
				const SequenceFactory &seqfac);
	
	void Process(Sequence &seq, double dTime);
	
private:
	IndelDistribution rand_ins;
	IndelDistribution rand_del;
	double m_dInsRate;
	double m_dDelRate;

	SequenceFactory rand_seq;
};

}

#endif
