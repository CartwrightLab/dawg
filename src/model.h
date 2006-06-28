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

	enum SeqType { TypeDNA, TypeRNA, TypeProtein, TypeCodon };

protected:
	bool SetupSubstModel(const std::string &ssModel, std::vector<double> &vdFreqs, std::vector<double> &vdParams, SeqType &tySeq);

	SubstModel m_modSubst;
	GillespieProcessor m_modIndel;
};

} // namespace dawg

#endif //DAWG_MODEL_H
