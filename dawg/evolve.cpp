#include "dawg.h"
#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <functional>
#include "rand.h"
#include "matrix.h"

#ifdef _WIN32
#	define finite _finite
#endif

using namespace std;

//void PrintMatrix(const SPMatrix& mat);

double g_dGamma = 0.0;  // The parameter of the Gamma distribution with mean 1
double g_dIota = 0.0;   // The probability that a site is invariant

double g_dNucCumFreqs[4] = {0.25, 0.5, 0.75, 1.0}; // ACGT cummulative frequencies.
double g_dFreqs[4] = {0.25, 0.25, 0.25, 0.25}; // ACGT Frequencies
Matrix44 g_matV, g_matU, g_matQ, g_matR;  // U = transpose(V)
Vector4  g_vecL;

/*string EvoDescription()
{
	ostringstream ss;
	ss << "Substitution Parameters:" << endl;
	ss << "  GTR Rate Matrix:" << endl;
	for(int i=0;i<4;++i)
	{
		ss << "    " << NucToChar((Nucleotide::Nuc)i) << " [";
		for(int j=0;j<4;++j)
			ss << g_matR(i,j) << ((j!=3) ? " " : "");
		ss << ']' << endl;
	}
	ss << "  Nucleotide Frequencies:" << endl;
	for(int i=0;i<4;++i)
		ss << "    " << NucToChar((Nucleotide::Nuc)i) << " = " << g_dFreqs[i] << endl;
	ss << "  Heterogenous Rates:" << endl;
	ss << "    Gamma = " << g_dGamma << endl;
	ss << "    Iota  = " << g_dIota << endl;
	ss << endl << "Gap Parameters:" << endl;
	ss << "  Insertions:" << endl;
	ss << "    Lambda = " << ((g_dLambda > DBL_EPSILON) ? g_sIndelParams.sIns.dLambda*g_dLambda : 0.0) << endl;
	if(g_bNegBnIndel)
		ss << "    Model = Negative Binomial(" << g_sIndelParams.sIns.uR << ", " << g_sIndelParams.sIns.dQ << ")" << endl;
	else
	{
		ss << "    Model = User(";
		vector<double>::const_iterator cit = g_sIndelParams.sIns.vSizes.begin();
		double d = *cit++;
		ss << d;
		for(;cit != g_sIndelParams.sIns.vSizes.end();++cit)
		{
			ss << ", " << *cit-d;
			d = *cit;
		}
		ss << ")" << endl;
	}
	ss << "  Deletions:" << endl;
	ss << "    Lambda = " << ((g_dLambda > DBL_EPSILON) ? g_sIndelParams.sDel.dLambda*g_dLambda : 0.0) << endl;
	if(g_bNegBnIndel)
		ss << "    Model = Negative Binomial( " << g_sIndelParams.sDel.uR << " " << g_sIndelParams.sDel.dQ << " )" << endl;
	else
	{
		ss << "    Model = User(";
		vector<double>::const_iterator cit = g_sIndelParams.sDel.vSizes.begin();
		double d = *cit++;
		ss << d;
		for(;cit != g_sIndelParams.sDel.vSizes.end();++cit)
		{
			ss << ", " << *cit-d;
			d = *cit;
		}
		ss << ")" << endl;
	}
	ss << endl << "Other Parameters:" << endl;
	ss << "   Branch Scale = " << g_dScale << endl;
	return ss.str();
}
*/

// Substitution Model is REV	
bool EvoRevParams(double pFreqs[], double pSubs[])
{
	if(pFreqs[0] < 0 || pFreqs[1] < 0 || pFreqs[2] < 0 || pFreqs[3] < 0)
		return DawgError("Nucleotide frequences need to be positive.");
	memcpy(g_dFreqs, pFreqs, 4*sizeof(double));
	g_dNucCumFreqs[0] = pFreqs[0];
	g_dNucCumFreqs[1] = pFreqs[1]+g_dNucCumFreqs[0];
	g_dNucCumFreqs[2] = pFreqs[2]+g_dNucCumFreqs[1];
	g_dNucCumFreqs[3] = 1.0;
	if(g_dNucCumFreqs[0] > 1.0 || g_dNucCumFreqs[1] > 1.0 || g_dNucCumFreqs[2] > 1.0)
		return DawgError("Nucleotide frequences need to sum to 1.");
	
	if(pSubs[0] < 0.0 || pSubs[1] < 0.0 || pSubs[2] < 0.0
		|| pSubs[3] < 0.0 || pSubs[4] < 0.0 || pSubs[5] < 0.0)
		return DawgError("Substitution rates need to be positive.");

	// Setup Symetric Matrix
	Matrix44 matQ(Matrix44::s_Zero);
	matQ(0,1) = matQ(1,0) = pSubs[0]; //A-C
	matQ(0,2) = matQ(2,0) = pSubs[1]; //A-G
	matQ(0,3) = matQ(3,0) = pSubs[2]; //A-T
	matQ(1,2) = matQ(2,1) = pSubs[3]; //C-G
	matQ(1,3) = matQ(3,1) = pSubs[4]; //C-T
	matQ(2,3) = matQ(3,2) = pSubs[5]; //G-T
	
	// Store Rate Matrix
	g_matR = matQ;
	
	// Create GTR Genetating Matrix
	Vector4 vecF(pFreqs);	
	matQ.Scale(matQ, vecF);

	// Scale such that the total rate of substitution is equal to one
	matQ(0,0) = -(matQ(0,1)+matQ(0,2)+matQ(0,3));
	matQ(1,1) = -(matQ(1,0)+matQ(1,2)+matQ(1,3));
	matQ(2,2) = -(matQ(2,0)+matQ(2,1)+matQ(2,3));
	matQ(3,3) = -(matQ(3,0)+matQ(3,1)+matQ(3,2));
	matQ.Scale(matQ, -1.0/(vecF[0]*matQ(0,0)+vecF[1]*matQ(1,1)+
		vecF[2]*matQ(2,2)+vecF[3]*matQ(3,3)));
	
	// Store Scaled Q Matrix
	g_matQ = matQ;

	// Make Q a symetric matrix again
	Vector4 vecD, vecE;  //D*E=I
	for(Matrix44::Pos i=0;i<4;++i)
	{
		vecD[i] = sqrt(vecF[i]);
		vecE[i] = 1.0/vecD[i];
	}
	matQ.Scale(matQ, vecE);
	matQ.Scale(vecD, matQ);

	//Find EigenSystem using Jacobian Transformations
	int nRet = EigenSystem(matQ, g_vecL, g_matV);
	if(nRet == -1)
		return DawgError("Eigensystem failed to converge.");
	g_matU = g_matV;
	g_matV.Scale(vecE, g_matV);
	g_matU.Transpose();
	g_matU.Scale(g_matU, vecD);

	Matrix44 matTemp1, matTemp2;
	matTemp1.Scale(g_vecL, g_matU);
	matTemp2.Multiply(g_matV,matTemp1);
	return true;
}

// Rates vary based on Gamma+Inv model
bool EvoRateParams(double dGamma, double dIota)
{
	bool bRet = true;
	if(dGamma < 0.0)
		bRet = DawgError("Invalid Gamma, \"%f\", Gamma must be positive.", dGamma);
	else if(0.0 > dIota || dIota > 1.0)
		bRet = DawgError("Invalid Iota, \"%f\", Iota must be a probability.", dIota);
	g_dGamma = dGamma;
	g_dIota = dIota;
	return bRet;
}

bool EvoScaleTrees(double dScale)
{
	if(!finite(dScale))
		return DawgError("Scale parameter needs to be finite.");
	g_dScale = dScale;
	return true;
}

// Gap model
bool EvoIndelUser(double dInsL, double pIns[], int nIns,
				  double dDelL, double pDel[], int nDel)
{
/*	g_bNegBnIndel = false;
	g_dLambda = dInsL+dDelL;
	g_sIndelParams.sIns.dLambda = dInsL/g_dLambda;
	g_sIndelParams.sDel.dLambda = dDelL/g_dLambda;

	g_sIndelParams.sIns.vSizes.clear();
	double dSum = 0.0;
	for(int i=0;i<nIns;++i)
	{
		if(pIns[i] > 1.0 || dSum > 1.0)
			return DawgError("User Gap Model (Ins) parameters must sum to one.");
		dSum += pIns[i];
		g_sIndelParams.sIns.vSizes.push_back(dSum);
	}
	g_sIndelParams.sIns.vSizes.back() = 1.0;
	g_sIndelParams.sDel.vSizes.clear();
	dSum = 0.0;
	for(int i=0;i<nDel;++i)
	{
		if(pDel[i] > 1.0 || dSum > 1.0)
			return DawgError("User Gap Model (Del) parameters must sum to one.");
		dSum += pDel[i];
		g_sIndelParams.sDel.vSizes.push_back(dSum);
	}
	g_sIndelParams.sDel.vSizes.back() = 1.0;
	*/
	return true;
}

bool EvoIndelNegBn(double dInsL, unsigned long uInsR, double dInsQ,
				   double dDelL, unsigned long uDelR, double dDelQ)
{
/*	bool bRet = true;
	g_bNegBnIndel = true;
	if(0.0 > dInsL || dInsL > 1.0)
		bRet = DawgError("NB Gap Model (Ins) Lambda must be a probability.");
	if(0.0 > dDelL || dDelL > 1.0)
		bRet = DawgError("NB Gap Model (Del) Lambda must be a probability.");
	if(0.0 > dInsQ || dInsQ > 1.0)
		bRet = DawgError("NB Gap Model (Ins) Q must be a probability.");
	if(0.0 > dDelQ || dDelQ > 1.0)
		bRet = DawgError("NB Gap Model (Del) Q must be a probability.");
	if(uInsR == 0)
		bRet = DawgError("NB Gap Model (Ins) R must be a positive integer.");
	if(uDelR == 0)
		bRet = DawgError("NB Gap Model (Del) R must be a positive integer.");

	g_dLambda = dInsL+dDelL;
	g_sIndelParams.sIns.dLambda = dInsL/g_dLambda;
	g_sIndelParams.sIns.uR = uInsR;
	g_sIndelParams.sIns.dQ = dInsQ;

	g_sIndelParams.sDel.dLambda = dDelL/g_dLambda;
	g_sIndelParams.sDel.uR = uDelR;
	g_sIndelParams.sDel.dQ = dDelQ; 
	return bRet;*/
	return true;
}

