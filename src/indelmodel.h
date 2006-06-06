#ifndef DAWG_INDELMODEL_H
#define DAWG_INDELMODEL_H

#include "rand.h"

namespace Dawg {

class IndelModel
{
public:
	typedef unsigned int result_type;
	
	enum Type {modZipf, modGeo, modUser};

	IndelModel();

	bool Create(Type mod, const std::vector<double> &vdParams);

	result_type Rand();
	
private:
	Type _model;
	std::vector<double> _vdParams;

	//zipf_distribution<result_type> zipf_dist;
	truncated_zipf_distribution<result_type> zipf_dist;
	discrete_distribution<result_type> user_dist;
	boost::geometric_distribution<result_type> geo_dist;

	boost::variate_generator<DawgRng&, truncated_zipf_distribution<result_type> > rand_zipf;
	boost::variate_generator<DawgRng&, discrete_distribution<result_type> > rand_user;
	boost::variate_generator<DawgRng&, boost::geometric_distribution<result_type> > rand_geo;
};

}

#endif
