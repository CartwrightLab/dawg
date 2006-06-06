#include "model.h"
#include <memory.h>

using namespace std;
//using namespace Dawg;

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

