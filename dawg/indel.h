#ifndef DAWG_INDEL_H
#define DAWG_INDEL_H

#include "dawg.h"

class Node;

class IndelModel
{
public:
	class Params
	{
	public:
		std::string ssModel;
		double dLambda;
		std::vector<double> vdModel;
	};
	virtual unsigned long RandSize() const = 0;
	virtual double MeanSize() const = 0;
};

class NegBnModel : public IndelModel
{
public:
	NegBnModel(const std::vector<double>& vdModel);
	virtual unsigned long RandSize() const;
	virtual double MeanSize() const;

protected:
	unsigned long m_uR;
	double m_dQ;
};

class UserModel : public IndelModel
{
public:
	UserModel(const std::vector<double>& vdModel);
	virtual unsigned long RandSize() const;
	virtual double MeanSize() const;

protected:
	std::vector<double> m_vSizesCum;
};

class LinearFunc : public std::unary_function<double, double>
{
public:
	result_type operator()( argument_type x) { return m*x+b; }
	argument_type m;
	argument_type b;
};

class IndelProcessor
{
public:
	bool Setup(const IndelModel::Params& rIns, const IndelModel::Params& rDel);
	void Process(Node* pNode);

protected:	
	std::auto_ptr<IndelModel> m_pInsertionModel;
	std::auto_ptr<IndelModel> m_pDeletionModel;
	double m_dLambdaIns;
	double m_dLambdaDel;
	LinearFunc m_funcRateIns;
	LinearFunc m_funcRateSum;
};

#endif // DAWG_INDEL_H