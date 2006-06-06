#ifndef DAWG_SUBSTMODEL_H
#define DAWG_SUBSTMODEL_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "sequence.h"

//namespace ublas = boost::numeric::ublas;

namespace Dawg {

class SubstModel
{
public:
	typedef Dawg::Sequence::base_type base_type;
	SubstModel() : dOldTime(-1.0) { }

	bool Create(const std::vector<double> &vdFreqs, const std::vector<double> &vdParams);
	void UpdateSubMatrix(double dTime);

	base_type Subst(base_type uCurr);
	base_type Subst(base_type uCurr, double dTime)
	{
		UpdateSubMatrix(dTime);
		return Subst(uCurr);
	}
	base_type operator()(base_type uCurr, double dTime) { return Subst(uCurr, dTime); }
	base_type operator()(base_type uCurr) { return Subst(uCurr); }
	
protected:

	ublas::matrix<double> matP;
	double dOldTime;

	ublas::matrix<double> matU;
	ublas::matrix<double> matV;
	ublas::vector<double> vecE;
};

}; // namespace Dawg

#endif
