#ifndef DAWG_INDELMODEL_H
#define DAWG_INDELMODEL_H

namespace Dawg {

class IndelModel
{
public:
	IndelModel();

	bool Create(const std::string& ssModel, const std::vector<double> &vdParams);
};

}

#endif
