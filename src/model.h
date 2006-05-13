#ifndef DAWG_MODEL_H
#define DAWG_MODEL_H

#include "dawgvar.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

namespace ublas = boost::numeric::ublas;

namespace Dawg {

class Model
{
public:
	Model();
	virtual ~Model();

	bool Create(const Dawg::Variables& var);
	bool Run();

protected:
	void UpdateSubMatrix(double dTime);
	ublas::matrix<double> matP;
	double dOldTime;

	ublas::matrix<double> matGTR;
	ublas::matrix<double> matU;
	ublas::matrix<double> matV;
	ublas::vector<double> vecE;

	double dSubFreqs[4];
	double dSubFreqsCum[4];

	double dSubParams[6];
	std::vector<double> dSubIota;
	std::vector<double> dSubScale;
	std::vector<double> dSubGamma;
};

} // namespace dawg

#endif //DAWG_MODEL_H
