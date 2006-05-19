#include "model.h"
#include <memory.h>
#include <gsl/gsl_eigen.h>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
//using namespace Dawg;

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
	return (status != 0);
}

bool Dawg::SubstModel::Create(const std::string &model, const double &s[6], const double &f[4],
				const std::vector<double> &iota, const std::vector<double> & scale)
{
	return true;
}

bool Dawg::SubstModel::Create(const std::vector<double> &f, const std::vector<double> &s,
	const std::vector<double> &iota, const std::vector<double> &scale)
{
	std::vector<double>::size_type N = f.size();
	if(N*(N-1)/2 != f.size())
		return DawgError("Sub.Params must be of length N(N-1)/2 where N is the length of Sub.Freqs.");
	double dTemp = 0.0;
	for(std::vector<double>::const_iterator cit = f.begin(); cit != f.end(); ++cit)
	{
		if(*cit < 0.0)
			return DawgError("Sub.Freqs must contain positive numbers.");
		dTemp += 0.0;
	}
	dTemp -= 1.0;
	if(dTemp > DBL_EPSILON || dTemp < -DBL_EPSILON)
		return DawgError("Sub.Freqs must sum to 1.0.");
	for(std::vector<double>::const_iterator cit = s.begin(); cit != s.end(); ++cit)
	{
		if(*cit < 0.0)
			return DawgError("Sub.Params must contain positive numbers.");
	}
	if(iota.size() != scale.size())
		return DawgError("Sub.Iota and Sub.Scale must have the same size.");
	
	double dY = 0.0;
	ublas::matrix<double> matQ(ublas::zero_matrix<double>(N,N));
	size_t x = 0;
	for(ublas::matrix<double>::size_type i = 0; i < N-1; ++i)
	{
		for(ublas::matrix<double>::size_type j = i+1; j < N; ++j, ++x)
		{
			matQ(i,j) = f[j]*s[x];
			matQ(i,i) += -f[j]*s[x];
			matQ(j,i) = f[i]*s[x];
			matQ(j,j) += -f[i]*s[x];
		}
		dY += f[i]*matQ(i,i);
	}
	
	// Scale such that the total rate of substitution is equal to one
	double dX = 0.0;
	for(size_t i=0;i<iota.size();++i)
		dX += (1.0-iota[i])*scale[i];
	dX = -((double)iota.size())/(dX*dY);
	matQ /= dX;
	
	// Make Q a symetric matrix again
	for(size_t i=0;i<N;++i)
	{
		double d = sqrt(v(i));
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
		double d = sqrt(v(i));
		for(size_t j=0;j<4;++j)
		{
			matV(i,j) /= d;
			matU(j,i) *= d;
		}
	}
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

Dawg::SubstModel::index_type Dawg::SubstModel::Subst(index_type p, double dTime)
{
	UpdateSubMatrix(dTime);
	return p;
}

bool Dawg::Model::Create(const Dawg::Variables& var)
{
	// Verifiy Parameters
	memcpy(dSubFreqs, var.dSubFreqs, sizeof(dSubFreqs));
	memcpy(dSubParams, var.dSubParams, sizeof(dSubParams));
	
	// Rate Heterogenity 
	dSubIota = var.dSubIota;
	dSubScale = var.dSubScale;
	dSubGamma = var.dSubGamma;

	// Cummulative Frequencies
	dSubFreqsCum[0] = dSubFreqs[0];
	dSubFreqsCum[1] = dSubFreqs[1]+dSubFreqsCum[0];
	dSubFreqsCum[2] = dSubFreqs[2]+dSubFreqsCum[1];
	dSubFreqsCum[3] = 1.0;

	// Create GTR Genetating Matrix
	if(!modSubst.Create(var.ssSubModel, var.dSubParams, var.dSubFreqs, var.dSubIota, var.dSubScale))
		return DawgError("Creation of substitution model failed.");
	
	return true;
}

