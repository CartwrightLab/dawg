#include "model.h"
#include <memory.h>
#include <gsl/gsl_eigen.h>
#include <boost/numeric/ublas/symmetric.hpp>

using namespace std;
//using namespace Dawg;

bool EigenSymmv(ublas::matrix<double> &A, ublas::vector<double> &E, ublas::matrix<double> &V)
{
	gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t n = A.size1();
	E.resize(4, false);
	V.resize(4, 4, false);
	gsl_matrix_view A_view = gsl_matrix_view_array(&A(0,0), n, n);
	gsl_vector_view E_view = gsl_vector_view_array(&E(0), n);
	gsl_matrix_view V_view = gsl_matrix_view_array(&V(0,0), n, n);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(4);

	int status = gsl_eigen_symmv(&A_view.matrix, &E_view.vector, &V_view.matrix, w);
	//gsl_eigen_symmv_sort(&E_view.vector, &V_view.matrix, GSL_EIGEN_SORT_ABS_ASC);
	gsl_eigen_symmv_free(w);
	gsl_set_error_handler(handler);
	return (status != 0);
}

bool Dawg::Model::Create(const Dawg::Variables& var)
{
	// Verifiy Parameters
	memcpy(dSubFreqs, var.dSubFreqs, sizeof(dSubFreqs));
	memcpy(dSubParams, var.dSubParams, sizeof(dSubParams));
	if(dSubFreqs[0] < 0.0 || dSubFreqs[1] < 0.0 || dSubFreqs[2] < 0.0 || dSubFreqs[3] < 0.0)
		return DawgError("Nucleotide frequences need to be positive.");
	dSubFreqs[3] = 1.0-dSubFreqs[0]-dSubFreqs[1]-dSubFreqs[2];
	if( var.dSubFreqs[3] < 0.0 )
		return DawgError("Nucleotide frequencies need to sum to 1.0.");
	if(dSubParams[0] < 0.0 || dSubParams[1] < 0.0 || dSubParams[2] < 0.0
		|| dSubParams[3] < 0.0 || dSubParams[4] < 0.0 || dSubParams[5] < 0.0)
		return DawgError("Substitution rates need to be positive.");
	
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
	ublas::matrix<double> matQ(ublas::zero_matrix<double>(4,4));
	size_t x = 0;
	double dY = 0.0;
	for(ublas::matrix<double>::size_type i = 0; i < matQ.size1()-1; ++i)
	{
		for(ublas::matrix<double>::size_type j = i+1; j < matQ.size2(); ++j, ++x)
		{
			matQ(i,j) = dSubFreqs[j]*dSubParams[x];
			matQ(i,i) += -dSubFreqs[j]*dSubParams[x];
			matQ(j,i) = dSubFreqs[i]*dSubParams[x];
			matQ(j,j) += -dSubFreqs[i]*dSubParams[x];
		}
		dY += dSubFreqs[i]*matQ(i,i);
	}
	
	// Scale such that the total rate of substitution is equal to one
	double dX = 0.0;
	for(size_t i=0;i<dSubIota.size();++i)
		dX += (1.0-dSubIota[i])*dSubScale[i];
	dX = -((double)dSubIota.size())/dX/dY;
	matQ /= dX;
	

	// Make Q a symetric matrix again
	for(size_t i=0;i<4;++i)
	{
		double d = sqrt(dSubFreqs[i]);
		for(size_t j=0;j<4;++j)
		{
			matQ(i,j) *= d;
			matQ(j,i) /= d;
		}
	}
	//Find EigenSystem
	if(!EigenSymmv(matQ, vecE, matV))
		return DawgError("Failure to calculate eigensystem");

	// Reset unsymetry
	matU = ublas::trans(matV);
	for(size_t i=0;i<4;++i)
	{
		double d = sqrt(dSubFreqs[i]);
		for(size_t j=0;j<4;++j)
		{
			matV(i,j) /= d;
			matU(j,i) *= d;
		}
	}

	return true;
}

void Dawg::Model::UpdateSubMatrix(double dTime)
{
	if(dTime == dOldTime)
		return;
	ublas::matrix<double> matT(ublas::zero_matrix<double>(4,4));
	for(ublas::matrix<double>::size_type i=0;i<matT.size1();++i)
		matT(i,i) = exp(dTime*vecE(i));
	noalias(matP) = ublas::prod(matV, ublas::matrix<double>(ublas::prod(matT, matU)));
	for(ublas::matrix<double>::size_type i=0;i<matT.size1();++i)
	{
		for(ublas::matrix<double>::size_type j=1;j<matT.size2()-1;++j)
			matP(i,j) += matP(i,j-1);
		//matP(i,matT.size2()) = 1.0;
	}
}