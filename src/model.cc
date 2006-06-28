#pragma warning(disable: 4127 4512)

#include "model.h"
#include <memory.h>

using namespace std;
//using namespace Dawg;

template<class T, class A, size_t N>
const T& assign_from_array(T& t, const A (&a)[N])
{
	t.assign(&a[0], &a[N]);
	return t;
}

double paramsWAG[] = {
	/*A-C*/ 1.027040,
	/*A-D*/ 0.738998,
	/*A-E*/ 1.582850,
	/*A-F*/ 0.210494,
	/*A-G*/ 1.416720,
	/*A-H*/ 0.316954,
	/*A-I*/ 0.193335,
	/*A-K*/ 0.906265,
	/*A-L*/ 0.397915,
	/*A-M*/ 0.893496,
	/*A-N*/ 0.509848,
	/*A-P*/ 1.438550,
	/*A-Q*/ 0.908598,
	/*A-R*/ 0.551571,
	/*A-S*/ 3.370790,
	/*A-T*/ 2.121110,
	/*A-V*/ 2.006010,
	/*A-W*/ 0.113133,
	/*A-Y*/ 0.240735,
	/*C-D*/ 0.0302949,
	/*C-E*/ 0.021352,
	/*C-F*/ 0.398020,
	/*C-G*/ 0.306674,
	/*C-H*/ 0.248972,
	/*C-I*/ 0.170135,
	/*C-K*/ 0.0740339,
	/*C-L*/ 0.384287,
	/*C-M*/ 0.390482,
	/*C-N*/ 0.265256,
	/*C-P*/ 0.109404,
	/*C-Q*/ 0.0988179,
	/*C-R*/ 0.528191,
	/*C-S*/ 1.407660,
	/*C-T*/ 0.512984,
	/*C-V*/ 1.002140,
	/*C-W*/ 0.717070,
	/*C-Y*/ 0.543833,
	/*D-E*/ 6.174160,
	/*D-F*/ 0.0467304,
	/*D-G*/ 0.865584,
	/*D-H*/ 0.930676,
	/*D-I*/ 0.039437,
	/*D-K*/ 0.479855,
	/*D-L*/ 0.0848047,
	/*D-M*/ 0.103754,
	/*D-N*/ 5.429420,
	/*D-P*/ 0.423984,
	/*D-Q*/ 0.616783,
	/*D-R*/ 0.147304,
	/*D-S*/ 1.071760,
	/*D-T*/ 0.374866,
	/*D-V*/ 0.152335,
	/*D-W*/ 0.129767,
	/*D-Y*/ 0.325711,
	/*E-F*/ 0.0811339,
	/*E-G*/ 0.567717,
	/*E-H*/ 0.570025,
	/*E-I*/ 0.127395,
	/*E-K*/ 2.584430,
	/*E-L*/ 0.154263,
	/*E-M*/ 0.315124,
	/*E-N*/ 0.947198,
	/*E-P*/ 0.682355,
	/*E-Q*/ 5.469470,
	/*E-R*/ 0.439157,
	/*E-S*/ 0.704939,
	/*E-T*/ 0.822765,
	/*E-V*/ 0.588731,
	/*E-W*/ 0.156557,
	/*E-Y*/ 0.196303,
	/*F-G*/ 0.049931,
	/*F-H*/ 0.679371,
	/*F-I*/ 1.059470,
	/*F-K*/ 0.088836,
	/*F-L*/ 2.115170,
	/*F-M*/ 1.190630,
	/*F-N*/ 0.0961621,
	/*F-P*/ 0.161444,
	/*F-Q*/ 0.0999208,
	/*F-R*/ 0.102711,
	/*F-S*/ 0.545931,
	/*F-T*/ 0.171903,
	/*F-V*/ 0.649892,
	/*F-W*/ 1.529640,
	/*F-Y*/ 6.454280,
	/*G-H*/ 0.249410,
	/*G-I*/ 0.0304501,
	/*G-K*/ 0.373558,
	/*G-L*/ 0.0613037,
	/*G-M*/ 0.174100,
	/*G-N*/ 1.125560,
	/*G-P*/ 0.243570,
	/*G-Q*/ 0.330052,
	/*G-R*/ 0.584665,
	/*G-S*/ 1.341820,
	/*G-T*/ 0.225833,
	/*G-V*/ 0.187247,
	/*G-W*/ 0.336983,
	/*G-Y*/ 0.103604,
	/*H-I*/ 0.138190,
	/*H-K*/ 0.890432,
	/*H-L*/ 0.499462,
	/*H-M*/ 0.404141,
	/*H-N*/ 3.956290,
	/*H-P*/ 0.696198,
	/*H-Q*/ 4.294110,
	/*H-R*/ 2.137150,
	/*H-S*/ 0.740169,
	/*H-T*/ 0.473307,
	/*H-V*/ 0.118358,
	/*H-W*/ 0.262569,
	/*H-Y*/ 3.873440,
	/*I-K*/ 0.323832,
	/*I-L*/ 3.170970,
	/*I-M*/ 4.257460,
	/*I-N*/ 0.554236,
	/*I-P*/ 0.0999288,
	/*I-Q*/ 0.113917,
	/*I-R*/ 0.186979,
	/*I-S*/ 0.319440,
	/*I-T*/ 1.458160,
	/*I-V*/ 7.821300,
	/*I-W*/ 0.212483,
	/*I-Y*/ 0.420170,
	/*K-L*/ 0.257555,
	/*K-M*/ 0.934276,
	/*K-N*/ 3.012010,
	/*K-P*/ 0.556896,
	/*K-Q*/ 3.894900,
	/*K-R*/ 5.351420,
	/*K-S*/ 0.967130,
	/*K-T*/ 1.386980,
	/*K-V*/ 0.305434,
	/*K-W*/ 0.137505,
	/*K-Y*/ 0.133264,
	/*L-M*/ 4.854020,
	/*L-N*/ 0.131528,
	/*L-P*/ 0.415844,
	/*L-Q*/ 0.869489,
	/*L-R*/ 0.497671,
	/*L-S*/ 0.344739,
	/*L-T*/ 0.326622,
	/*L-V*/ 1.800340,
	/*L-W*/ 0.665309,
	/*L-Y*/ 0.398618,
	/*M-N*/ 0.198221,
	/*M-P*/ 0.171329,
	/*M-Q*/ 1.545260,
	/*M-R*/ 0.683162,
	/*M-S*/ 0.493905,
	/*M-T*/ 1.516120,
	/*M-V*/ 2.058450,
	/*M-W*/ 0.515706,
	/*M-Y*/ 0.428437,
	/*N-P*/ 0.195081,
	/*N-Q*/ 1.543640,
	/*N-R*/ 0.635346,
	/*N-S*/ 3.974230,
	/*N-T*/ 2.030060,
	/*N-V*/ 0.196246,
	/*N-W*/ 0.0719167,
	/*N-Y*/ 1.086000,
	/*P-Q*/ 0.933372,
	/*P-R*/ 0.679489,
	/*P-S*/ 1.613280,
	/*P-T*/ 0.795384,
	/*P-V*/ 0.314887,
	/*P-W*/ 0.139405,
	/*P-Y*/ 0.216046,
	/*Q-R*/ 3.035500,
	/*Q-S*/ 1.028870,
	/*Q-T*/ 0.857928,
	/*Q-V*/ 0.301281,
	/*Q-W*/ 0.215737,
	/*Q-Y*/ 0.227710,
	/*R-S*/ 1.224190,
	/*R-T*/ 0.554413,
	/*R-V*/ 0.251849,
	/*R-W*/ 1.163920,
	/*R-Y*/ 0.381533,
	/*S-T*/ 4.378020,
	/*S-V*/ 0.232739,
	/*S-W*/ 0.523742,
	/*S-Y*/ 0.786993,
	/*T-V*/ 1.388230,
	/*T-W*/ 0.110864,
	/*T-Y*/ 0.291148,
	/*V-W*/ 0.365369,
	/*V-Y*/ 0.314730,
	/*W-Y*/ 2.485390,
};

double freqsWAG[] = {
	/*A*/ 0.0866279,
	/*C*/ 0.0193078,
	/*D*/ 0.0570451,
	/*E*/ 0.0580589,
	/*F*/ 0.0384319,
	/*G*/ 0.0832518,
	/*H*/ 0.0244313,
	/*I*/ 0.048466,
	/*K*/ 0.0620286,
	/*L*/ 0.086209,
	/*M*/ 0.0195027,
	/*N*/ 0.0390894,
	/*P*/ 0.0457631,
	/*Q*/ 0.0367281,
	/*R*/ 0.043972,
	/*S*/ 0.0695179,
	/*T*/ 0.0610127,
	/*V*/ 0.0708956,
	/*W*/ 0.0143859,
	/*Y*/ 0.0352742,
};

double paramsWAGStar[] = {
	/*A-C*/ 1.21324,
	/*A-D*/ 0.731152,
	/*A-E*/ 1.55788,
	/*A-F*/ 0.213179,
	/*A-G*/ 1.41993,
	/*A-H*/ 0.317684,
	/*A-I*/ 0.214596,
	/*A-K*/ 0.881639,
	/*A-L*/ 0.400822,
	/*A-M*/ 0.887458,
	/*A-N*/ 0.514347,
	/*A-P*/ 1.51861,
	/*A-Q*/ 1.03344,
	/*A-R*/ 0.589718,
	/*A-S*/ 3.52499,
	/*A-T*/ 2.24161,
	/*A-V*/ 1.92496,
	/*A-W*/ 0.135395,
	/*A-Y*/ 0.270321,
	/*C-D*/ 0.0379056,
	/*C-E*/ 0.0284956,
	/*C-F*/ 0.485001,
	/*C-G*/ 0.312544,
	/*C-H*/ 0.341479,
	/*C-I*/ 0.198958,
	/*C-K*/ 0.0719929,
	/*C-L*/ 0.451124,
	/*C-M*/ 0.428648,
	/*C-N*/ 0.233527,
	/*C-P*/ 0.109081,
	/*C-Q*/ 0.0999068,
	/*C-R*/ 0.568449,
	/*C-S*/ 1.35221,
	/*C-T*/ 0.522957,
	/*C-V*/ 1.10899,
	/*C-W*/ 0.728065,
	/*C-Y*/ 0.481954,
	/*D-E*/ 6.04299,
	/*D-F*/ 0.0458258,
	/*D-G*/ 0.88357,
	/*D-H*/ 0.958529,
	/*D-I*/ 0.0390513,
	/*D-K*/ 0.480308,
	/*D-L*/ 0.0869637,
	/*D-M*/ 0.0992829,
	/*D-N*/ 5.30821,
	/*D-P*/ 0.444152,
	/*D-Q*/ 0.657364,
	/*D-R*/ 0.159054,
	/*D-S*/ 1.09965,
	/*D-T*/ 0.395176,
	/*D-V*/ 0.155419,
	/*D-W*/ 0.142159,
	/*D-Y*/ 0.326191,
	/*E-F*/ 0.0873936,
	/*E-G*/ 0.588609,
	/*E-H*/ 0.599188,
	/*E-I*/ 0.124553,
	/*E-K*/ 2.45392,
	/*E-L*/ 0.154936,
	/*E-M*/ 0.294481,
	/*E-N*/ 1.00122,
	/*E-P*/ 0.720567,
	/*E-Q*/ 5.6037,
	/*E-R*/ 0.443685,
	/*E-S*/ 0.822025,
	/*E-T*/ 0.889765,
	/*E-V*/ 0.588443,
	/*E-W*/ 0.176397,
	/*E-Y*/ 0.209621,
	/*F-G*/ 0.0552962,
	/*F-H*/ 0.631713,
	/*F-I*/ 1.06458,
	/*F-K*/ 0.0832422,
	/*F-L*/ 2.10414,
	/*F-M*/ 1.14516,
	/*F-N*/ 0.0848492,
	/*F-P*/ 0.165205,
	/*F-Q*/ 0.109241,
	/*F-R*/ 0.122792,
	/*F-S*/ 0.563999,
	/*F-T*/ 0.188237,
	/*F-V*/ 0.653015,
	/*F-W*/ 1.58681,
	/*F-Y*/ 6.49269,
	/*G-H*/ 0.279542,
	/*G-I*/ 0.0310522,
	/*G-K*/ 0.381514,
	/*G-L*/ 0.067443,
	/*G-M*/ 0.184545,
	/*G-N*/ 1.12717,
	/*G-P*/ 0.254626,
	/*G-Q*/ 0.346823,
	/*G-R*/ 0.629768,
	/*G-S*/ 1.33618,
	/*G-T*/ 0.236489,
	/*G-V*/ 0.190095,
	/*G-W*/ 0.366467,
	/*G-Y*/ 0.108982,
	/*H-I*/ 0.162975,
	/*H-K*/ 0.854485,
	/*H-L*/ 0.508952,
	/*H-M*/ 0.40117,
	/*H-N*/ 3.9337,
	/*H-P*/ 0.722123,
	/*H-Q*/ 4.87366,
	/*H-R*/ 2.31211,
	/*H-S*/ 0.876688,
	/*H-T*/ 0.54992,
	/*H-V*/ 0.119749,
	/*H-W*/ 0.261223,
	/*H-Y*/ 4.31772,
	/*I-K*/ 0.320597,
	/*I-L*/ 3.1554,
	/*I-M*/ 3.94646,
	/*I-N*/ 0.527321,
	/*I-P*/ 0.111722,
	/*I-Q*/ 0.125999,
	/*I-R*/ 0.187262,
	/*I-S*/ 0.321774,
	/*I-T*/ 1.48876,
	/*I-V*/ 7.48376,
	/*I-W*/ 0.259584,
	/*I-Y*/ 0.44009,
	/*K-L*/ 0.255092,
	/*K-M*/ 0.877057,
	/*K-N*/ 2.88102,
	/*K-P*/ 0.588203,
	/*K-Q*/ 4.19125,
	/*K-R*/ 5.74119,
	/*K-S*/ 1.05314,
	/*K-T*/ 1.45173,
	/*K-V*/ 0.300343,
	/*K-W*/ 0.159261,
	/*K-Y*/ 0.155623,
	/*L-M*/ 4.81956,
	/*L-N*/ 0.144354,
	/*L-P*/ 0.422851,
	/*L-Q*/ 0.873266,
	/*L-R*/ 0.51821,
	/*L-S*/ 0.351913,
	/*L-T*/ 0.351564,
	/*L-V*/ 1.82105,
	/*L-W*/ 0.706082,
	/*L-Y*/ 0.427718,
	/*M-N*/ 0.198404,
	/*M-P*/ 0.179858,
	/*M-Q*/ 1.64018,
	/*M-R*/ 0.660816,
	/*M-S*/ 0.554077,
	/*M-T*/ 1.56873,
	/*M-V*/ 2.03324,
	/*M-W*/ 0.565299,
	/*M-Y*/ 0.437069,
	/*N-P*/ 0.204905,
	/*N-Q*/ 1.62299,
	/*N-R*/ 0.67416,
	/*N-S*/ 3.90127,
	/*N-T*/ 2.06787,
	/*N-V*/ 0.193323,
	/*N-W*/ 0.0746093,
	/*N-Y*/ 1.05269,
	/*P-Q*/ 0.913179,
	/*P-R*/ 0.711498,
	/*P-S*/ 1.54694,
	/*P-T*/ 0.802531,
	/*P-V*/ 0.325745,
	/*P-W*/ 0.135024,
	/*P-Y*/ 0.212945,
	/*Q-R*/ 3.02808,
	/*Q-S*/ 0.87908,
	/*Q-T*/ 0.829315,
	/*Q-V*/ 0.32893,
	/*Q-W*/ 0.208163,
	/*Q-Y*/ 0.210494,
	/*R-S*/ 1.35611,
	/*R-T*/ 0.594177,
	/*R-V*/ 0.282892,
	/*R-W*/ 1.24086,
	/*R-Y*/ 0.386714,
	/*S-T*/ 4.02507,
	/*S-V*/ 0.23769,
	/*S-W*/ 0.528249,
	/*S-Y*/ 0.742154,
	/*T-V*/ 1.4088,
	/*T-W*/ 0.118584,
	/*T-Y*/ 0.286443,
	/*V-W*/ 0.396884,
	/*V-Y*/ 0.353358,
	/*W-Y*/ 2.42261,
};

double freqsWAGStar[] = {
	/*A*/ 0.0866279,
	/*C*/ 0.0193078,
	/*D*/ 0.0570451,
	/*E*/ 0.0580589,
	/*F*/ 0.0384319,
	/*G*/ 0.0832518,
	/*H*/ 0.0244313,
	/*I*/ 0.048466,
	/*K*/ 0.0620286,
	/*L*/ 0.086209,
	/*M*/ 0.0195027,
	/*N*/ 0.0390894,
	/*P*/ 0.0457631,
	/*Q*/ 0.0367281,
	/*R*/ 0.043972,
	/*S*/ 0.0695179,
	/*T*/ 0.0610127,
	/*V*/ 0.0708956,
	/*W*/ 0.0143859,
	/*Y*/ 0.0352742,
};


bool Dawg::Model::SetupSubstModel(const std::string &ssModel, std::vector<double> &vdFreqs, std::vector<double> &vdParams, SeqType &tySeq)
{
	if(ssModel == "JC")
	{
		vdFreqs.assign(4, 0.25);
		vdParams.assign(6, 1.0);
		tySeq = TypeDNA;
	}
	else if(ssModel == "GTR")
	{
		if(vdFreqs.size() != 4 || vdParams.size() != 6)
			return DawgError("GTR model requires 4 frequencies and 6 numerical parameters.");
		tySeq = TypeDNA;
	}
	else if(ssModel == "K2P")
	{
		if(vdParams.size() < 2)
			return DawgError("K2P model requires 2 numerical parameters.");
		vdFreqs.assign(4, 0.25);
		double dTi = vdParams[0];
		double dTv = vdParams[1];
		vdParams.assign(6, dTv);
		vdParams[4] = vdParams[1] = dTi;
		tySeq = TypeDNA;
	}
	else if(ssModel == "K3P")
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
		tySeq = TypeDNA;
	}
	else if(ssModel == "HKY")
	{
		if(vdFreqs.size() != 4 || vdParams.size() < 1)
			return DawgError("HKY model requires 4 frequencies and 1 numerical parameter.");
		double dTv = vdParams[0];
		vdParams.assign(6, 1.0);
		vdParams[4] = vdParams[1] = dTv;
		tySeq = TypeDNA;
	}
	else if(ssModel == "F81")
	{
		if(vdFreqs.size() != 4)
			return DawgError("F81 model requires 4 frequencies.");
		vdParams.assign(6, 1.0);
		tySeq = TypeDNA;
	}
	else if(ssModel == "F84")
	{
		if(vdFreqs.size() != 4 || vdParams.size() < 1)
			return DawgError("F84 model requires 4 frequencies and 1 numerical parameter.");
		double dTv1 = 1.0 + vdParams[0]/(vdFreqs[0]+vdFreqs[2]);
		double dTv4 = 1.0 + vdParams[0]/(vdFreqs[1]+vdFreqs[3]);
		vdParams.assign(6, 1.0);
		vdParams[1] = dTv1;
		vdParams[4] = dTv4;
		tySeq = TypeDNA;
	}
	else if(ssModel == "TN")
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
		tySeq = TypeDNA;
	}
	else if(ssModel == "WAG")
	{
		assign_from_array(vdFreqs, freqsWAG);
		assign_from_array(vdParams, paramsWAG);
		tySeq = TypeProtein;
	}
	else if(ssModel == "WAG*")
	{
		assign_from_array(vdFreqs, freqsWAGStar);
		assign_from_array(vdParams, paramsWAGStar);
		tySeq = TypeProtein;
	}
	else
	{
		return DawgError("Unknown substitution model, \"%s\".", ssModel.c_str());
	}	
	return true;
}

bool Dawg::Model::Create(const Dawg::Variables& var)
{
	// Create Substitution Model
	std::vector<double> vdFreqs;
	std::vector<double> vdParams;
	SeqType typeSeq;
	SetupSubstModel(var.ssSubModel, vdFreqs, vdParams, typeSeq);
	if(!m_modSubst.Create(vdFreqs, vdParams))
		return DawgError("Creation of substitution model failed.");
	
	return true;
}

