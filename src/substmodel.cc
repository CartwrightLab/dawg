#include "dawg.h"
#include "substmodel.h"
#include "rand.h"
#include <gsl/gsl_eigen.h>

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

//bool Dawg::SubstModel::Create(const std::string &model, const std::vector<double> &f, const std::vector<double> &s,
//				const std::vector<double> &iota, const std::vector<double> & scale)
//{
//	std::vector<double> freqs(4, 0.25);	
//	std::vector<double> rates(6, 1.0);
//	if(model == "JC")
//	{
//		// do nothing
//	}
//	else if(model == "GTR")
//	{
//		if(f.size() < 4 || s.size() < 6)
//			return DawgError("GTR model requires 4 frequencies and 6 numerical parameters.");
//		freqs = f;
//		rates = s;
//	}
//	else if(model == "K2P")
//	{
//		if(s.size() < 2)
//			return DawgError("K2P model requires 2 numerical parameters.");
//		rates[4] = rates[1] = s[0];
//		rates[0] = rates[2] = rates[3] = rates[5] = s[1];
//	}
//	else if(model == "K3P")
//	{
//		if(s.size() < 2)
//			return DawgError("K2P model requires 3 numerical parameters.");
//		rates[4] = rates[1] = s[0];
//		rates[2] = rates[3] = s[1];
//		rates[0] = rates[5] = s[2];
//
//	}
//	else if(model == "HKY")
//	{
//		if(f.size() < 4 || s.size() < 1)
//			return DawgError("HKY model requires 4 frequencies and 1 numerical parameter.");
//		freqs = f;
//		rates[4] = rates[1] = s[0];
//		rates[0] = rates[2] = rates[3] = rates[5] = 1.0;
//
//	}
//	else if(model == "F81")
//	{
//		if(f.size() < 4)
//			return DawgError("F81 model requires 4 frequencies.");
//		freqs = f;
//	}
//	else if(model == "F84")
//	{
//		if(f.size() < 4 && s.size() < 1)
//			return DawgError("F84 model requires 4 frequencies and 1 numerical parameter.");
//		freqs = f;
//		rates[1] = 1.0 + s[0]/(f[0]+f[2]);
//		rates[4] = 1.0 + s[0]/(f[1]+f[3]);
//		rates[0] = rates[2] = rates[3] = rates[5] = 1.0;
//
//	}
//	else if(model == "TN")
//	{
//		if(vdParams.size() < 3)
//			return DawgError("TN model requires 4 frequencies and 3 numerical parameters.");
//		freqs = f;
//		rates[1] = s[0];
//		rates[4] = s[1];
//		rates[0] = rates[2] = rates[3] = rates[5] = s[2];
//	}
//	else
//	{
//		return DawgError("Unknown substitution model, \"%s\".", model.c_str());
//	}
//	return Create(freqs, rates, iota, scale);	
//}

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
	for(ublas::matrix<double>::size_type i = 0; i < N-1; ++i)
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
