#ifndef DAWG_MODEL_H
#define DAWG_MODEL_H

#include "dawgvar.h"
#include "sequence.h"
#include "substmodel.h"
#include "indelmodel.h"

namespace Dawg {

class Model
{
public:
	Model();
	virtual ~Model();

	bool Create(const Dawg::Variables& var);
	bool Run();

protected:
	bool SetupSubstModel(const std::string &ssModel, std::vector<double> &vdFreqs, std::vector<double> &vdParams, SeqType &tySeq);
	bool SetupInsModel(const std::string &ssModel, double dLambda, const std::vector<double> &vdParams);
	bool SetupDelModel(const std::string &ssModel, double dLambda, const std::vector<double> &vdParams);

	SubstModel m_modSubst;
	IndelModel m_modIndel;
	SequenceFactory m_facSeq;
};

} // namespace dawg

#endif //DAWG_MODEL_H
