#pragma warning(disable: 4127 4512)

#include "model.h"
#include <memory.h>

// Include Model Arrays
#include "modles.inl"

using namespace std;
//using namespace Dawg;



bool Dawg::Model::Create(const Dawg::Variables& var)
{
	// Create Substitution Model
	std::vector<double> vdFreqs;
	std::vector<double> vdParams;
	SeqType typeSeq;
	SetupSubstModel(var.ssSubModel, vdFreqs, vdParams, typeSeq);
	if(!m_modSubst.Create(vdFreqs, vdParams))
		return DawgError("Creation of substitution model failed.");
	if(!m_facSeq.Create(vdFreqs, var.dSubGamma, var.dSubIota))
		return DawgError("Creation of sequence model failed.");

	return true;
}

