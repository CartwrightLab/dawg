#include "dawg.h"
#include "dawgvar.h"
#include "substmodel.h"
#include "rand.h"
#include <gsl/gsl_eigen.h>

template<class T, class A, size_t N>
const T& assign_from_array(T& t, const A (&a)[N])
{
	t.assign(&a[0], &a[N]);
	return t;
}

bool EigenSymmv(ublas::matrix<double> &A, ublas::vector<double> &E, ublas::matrix<double> &V)
{
	gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t N = A.size1();
	E.resize(N, false);
	V.resize(N, N, false);
	gsl_matrix_view A_view = gsl_matrix_view_array(&A(0,0), N, N);
	gsl_vector_view E_view = gsl_vector_view_array(&E(0), N);
	gsl_matrix_view V_view = gsl_matrix_view_array(&V(0,0), N, N);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);

	int status = gsl_eigen_symmv(&A_view.matrix, &E_view.vector, &V_view.matrix, w);
	//gsl_eigen_symmv_sort(&E_view.vector, &V_view.matrix, GSL_EIGEN_SORT_ABS_ASC);
	gsl_eigen_symmv_free(w);
	gsl_set_error_handler(handler);
	return (status == 0);
}

bool Dawg::SubstModel::Create(const Dawg::Variables& var)
{
	vector<double> vdFreqs(var.dSubFreqs);
	vector<double> vdParams(var.dSubParams);

	if(var.ssSubModel == "JC")
	{
		vdFreqs.assign(4, 0.25);
		vdParams.assign(6, 1.0);
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "GTR")
	{
		if(vdFreqs.size() != 4 || vdParams.size() != 6)
			return DawgError("GTR model requires 4 frequencies and 6 numerical parameters.");
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "K2P")
	{
		if(vdParams.size() < 2)
			return DawgError("K2P model requires 2 numerical parameters.");
		vdFreqs.assign(4, 0.25);
		double dTi = vdParams[0];
		double dTv = vdParams[1];
		vdParams.assign(6, dTv);
		vdParams[4] = vdParams[1] = dTi;
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "K3P")
	{
		if(vdParams.size() < 2)
			return DawgError("K2P model requires 3 numerical parameters.");
		vdFreqs.assign(4, 0.25);
		double dTi = vdParams[0];
		double dTv1 = vdParams[1];
		double dTv2 = vdParams[2];
		vdParams.assign(6, dTi);
		vdParams[2] = vdParams[3] = dTv1;
		vdParams[0] = vdParams[5] = dTv2;
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "HKY")
	{
		if(vdFreqs.size() != 4 || vdParams.size() < 1)
			return DawgError("HKY model requires 4 frequencies and 1 numerical parameter.");
		double dTv = vdParams[0];
		vdParams.assign(6, 1.0);
		vdParams[4] = vdParams[1] = dTv;
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "F81")
	{
		if(vdFreqs.size() != 4)
			return DawgError("F81 model requires 4 frequencies.");
		vdParams.assign(6, 1.0);
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "F84")
	{
		if(vdFreqs.size() != 4 || vdParams.size() < 1)
			return DawgError("F84 model requires 4 frequencies and 1 numerical parameter.");
		double dTv1 = 1.0 + vdParams[0]/(vdFreqs[0]+vdFreqs[2]);
		double dTv4 = 1.0 + vdParams[0]/(vdFreqs[1]+vdFreqs[3]);
		vdParams.assign(6, 1.0);
		vdParams[1] = dTv1;
		vdParams[4] = dTv4;
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "TN")
	{
		if(vdParams.size() < 3)
			return DawgError("TN model requires 4 frequencies and 3 numerical parameters.");
		vdFreqs.assign(4, 0.25);
		double dTv1 = vdParams[0];
		double dTv4 = vdParams[1];
		double dTi = vdParams[2];
		vdParams.assign(6, dTi);
		vdParams[1] = dTv1;
		vdParams[4] = dTv4;
		m_type = TypeDNA;
	}
	else if(var.ssSubModel == "WAG")
	{
		assign_from_array(vdFreqs, freqsWAG);
		assign_from_array(vdParams, paramsWAG);
		m_type = TypeProtein;
	}
	else if(var.ssSubModel == "WAG*")
	{
		assign_from_array(vdFreqs, freqsWAGStar);
		assign_from_array(vdParams, paramsWAGStar);
		m_type = TypeProtein;
	}
	else
	{
		return DawgError("Unknown substitution model, \"%s\".", ssModel.c_str());
	}	
	return Create(vdFreqs, vdParams);
}

bool Dawg::SubstModel::Create(const std::vector<double> &vdFreqs, const std::vector<double> &vdParams)
{
	std::vector<double>::size_type N = vdFreqs.size();
	if(N*(N-1)/2 != vdParams.size())
		return DawgError("Sub.Params must be of length N(N-1)/2 where N is the length of Sub.Freqs.");
	double dTemp = 0.0;
	for(std::vector<double>::const_iterator cit = vdFreqs.begin(); cit != vdFreqs.end(); ++cit)
	{
		if(*cit < 0.0)
			return DawgError("Sub.Freqs must contain positive numbers.");
		dTemp += 0.0;
	}
	dTemp -= 1.0;
	if(dTemp > DBL_EPSILON || dTemp < -DBL_EPSILON)
		return DawgError("Sub.Freqs must sum to 1.0.");
	for(std::vector<double>::const_iterator cit = vdParams.begin(); cit != vdParams.end(); ++cit)
	{
		if(*cit < 0.0)
			return DawgError("Sub.Params must contain positive numbers.");
	}

	double dY = 0.0;
	ublas::matrix<double> matQ(ublas::zero_matrix<double>(N,N));
	size_t x = 0;
	for(ublas::matrix<double>::size_type i = 0; i < N; ++i)
	{
		for(ublas::matrix<double>::size_type j = i+1; j < N; ++j, ++x)
		{
			matQ(i,j) = vdFreqs[j]*vdParams[x];
			matQ(i,i) += -vdFreqs[j]*vdParams[x];
			matQ(j,i) = vdFreqs[i]*vdParams[x];
			matQ(j,j) += -vdFreqs[i]*vdParams[x];
		}
		dY += vdFreqs[i]*matQ(i,i);
	}
	
	// Scale such that the total rate of substitution is equal to one
	matQ /= dY;
	
	// Make Q a symetric matrix again
	for(size_t i=0;i<N;++i)
	{
		double d = sqrt(vdFreqs[i]);
		for(size_t j=0;j<N;++j)
		{
			matQ(i,j) *= d;
			matQ(j,i) /= d;
		}
	}
	//Find EigenSystem
	if(!EigenSymmv(matQ, vecE, matV))
		return DawgError("Failure to calculate eigensystem for substitution model.");

	// Reset unsymetry
	matU = ublas::trans(matV);
	for(size_t i=0;i<4;++i)
	{
		double d = sqrt(vdFreqs[i]);
		for(size_t j=0;j<4;++j)
		{
			matV(i,j) /= d;
			matU(j,i) *= d;
		}
	}
	dOldTime = -1.0;
	UpdateSubMatrix(1.0);
	return true;
}

void Dawg::SubstModel::UpdateSubMatrix(double dTime)
{
	if(dTime == dOldTime)
		return;
	ublas::matrix<double>::size_type N = vecE.size();
	ublas::matrix<double> matT(ublas::zero_matrix<double>(N,N));
	for(ublas::matrix<double>::size_type i=0;i<N;++i)
		matT(i,i) = exp(dTime*vecE(i));
	noalias(matP) = ublas::prod(matV, ublas::matrix<double>(ublas::prod(matT, matU)));
	for(ublas::matrix<double>::size_type i=0;i<N;++i)
	{
		for(ublas::matrix<double>::size_type j=1;j<N-1;++j)
			matP(i,j) += matP(i,j-1);
		//matP(i,matT.size2()) = 1.0;
	}
}

Dawg::SubstModel::base_type Dawg::SubstModel::Subst(base_type uCurr)
{
	double d = g_randReal01();
	const base_type N = vecE.size();
	base_type u;
	for(u=0; u < N-1 && d >= matP(uCurr, u); ++u) { }
	return u;
}

void Dawg::SubstModel::Process(Dawg::Sequence &seq, double dTime)
{
	for(Dawg::Sequence::pos_type u = 0; u < seq.Length(); ++u)
		seq.Residue(u).first = Subst(seq.Base(u), dTime*seq.Rate(u));
}

