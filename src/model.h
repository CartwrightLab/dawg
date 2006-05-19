#ifndef DAWG_MODEL_H
#define DAWG_MODEL_H

#include "dawgvar.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

namespace ublas = boost::numeric::ublas;

namespace Dawg {

class SubstModel
{
public:
	typedef size_t index_type;
	SubstModel() : dOldTime(1.0) { }
	void UpdateSubMatrix(double dTime);
	
	bool Create(const std::string &model,
				const double &s[6],
				const double &f[4],
				const std::vector<double> &iota,
				const std::vector<double> & scale);
	
	bool Create(const std::vector<double> &f, const std::vector<double> &s,
				const std::vector<double> &iota, const std::vector<double> &scale);
				
	index_type Subst(index_type p, double dTime);
	
protected:
	ublas::matrix<double> matP;
	double dOldTime;

	ublas::matrix<double> matU;
	ublas::matrix<double> matV;
	ublas::vector<double> vecE;
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
