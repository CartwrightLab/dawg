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

string EvoDescription()
{
	return "";
}

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

