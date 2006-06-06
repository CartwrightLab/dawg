#ifndef DAWG_MODEL_H
#define DAWG_MODEL_H

#include "dawgvar.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

namespace ublas = boost::numeric::ublas;

namespace Dawg {

class IndelModel
{
public:
	IndelModel();

	bool Create(const std::string& ssModel, const std::vector<double> &vdParams);
};

class Model
{
public:
	Model();
	virtual ~Model();

	bool Create(const Dawg::Variables& var);
	bool Run();

protected:
	SubstModel modSubst;
	
	double dSubFreqs[4];
	double dSubFreqsCum[4];

	double dSubParams[6];
	std::vector<double> dSubIota;
	std::vector<double> dSubScale;
	std::vector<double> dSubGamma;

};

} // namespace dawg

#endif //DAWG_MODEL_H
