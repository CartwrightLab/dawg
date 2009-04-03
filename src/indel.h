// indel.h - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_INDEL_H
#define DAWG_INDEL_H

#include "dawg.h"

class Node;

// Base class for indel models
class IndelModel
{
public:
	typedef std::vector<double>::size_type size_type;
	class Params
	{
	public:
		std::string ssModel;
		double dLambda;
		std::vector<double> vdModel;
	};
	virtual size_type RandSize() const = 0;
	virtual double MeanSize() const = 0;
	virtual ~IndelModel() { }
};

// Negative Binomial Model
// f(x) = NCh(r+x-2,x-1) (1-q)^r q^(x-1)
class NegBnModel : public IndelModel
{
public:
	NegBnModel(const std::vector<double>& vdModel);
	virtual size_type RandSize() const;
	virtual double MeanSize() const;

protected:
	unsigned int m_uR;
	double m_dQ;
};


// User Model
// f(x) = parameter(x-1)
class UserModel : public IndelModel
{
public:
	UserModel(const std::vector<double>& vdModel);
	virtual size_type RandSize() const;
	virtual double MeanSize() const;

protected:
	std::vector<double> m_vSizesCum;
	double m_dMean;
};

// Power-Law Model
// f(x) is proportional to (x^-a)
class PowerModel : public IndelModel
{
public:
	PowerModel(const std::vector<double>& vdModel);
	virtual size_type RandSize() const;
	virtual double MeanSize() const;
protected:
	double m_dA;
	size_type m_uM;
	double m_dMean;
};


// A class wrapping y = m*x + b
class LinearFunc : public std::unary_function<double, double>
{
public:
	result_type operator()( argument_type x) { return m*x+b; }
	argument_type m;
	argument_type b;
};

#endif // DAWG_INDEL_H

