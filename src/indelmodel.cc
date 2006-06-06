#include "indelmodel.h"

bool Dawg::IndelModel::Create(Type mod, const std::vector<double> &vdParams)
{
	switch(_model)
	{
	case modUser:
		user_dist = discrete_distribution<result_type>(vdParams.begin(), vdParams.end());
		return true;
	case modGeo:
		geo_dist = boost::geometric_distribution<result_type>(vdParams[0]);
		return true;
	case modZipf:
		zipf_dist = truncated_zipf_distribution<result_type>(vdParams[0], result_type(vdParams[1]));
		return true;
	default:
		return false;
	};
}

Dawg::IndelModel::result_type Dawg::IndelModel::Rand()
{
	switch(_model)
	{
	case modUser:
		return rand_user();
	case modGeo:
		return rand_geo();
	case modZipf:
	default:
		return rand_zipf();
	};
}

