#ifndef DAWG_INDELMODEL_H
#define DAWG_INDELMODEL_H

#include "rand.h"
#include "sequence.h"

#include <boost/ptr_container/ptr_vector.hpp>

namespace Dawg {

class IndelModel
{
public:
	class Element
	{
	public:
		virtual double Rate(const Sequence& /*seq*/) const { return 0.0;}
		virtual void operator()(Sequence& /*seq*/) { }
		void Process(Sequence & seq) { (*this)(seq); }
	};
	
	bool Create(const Dawg::Variables& var);
	bool Create(Element* pIns, Element* pDel)
	{
		AddElement(pIns);
		AddElement(pDel);
		return true;
	}

	void operator()(Sequence& seq, double dTime);

private:
	void AddElement(Element* p) { m_elements.push_back(p); }

	double Waiting(double d) { return rand_exp(d); }
	double Which(double d) { return rand_real(0.0, d); }

	boost::ptr_vector<Element> m_elements;
};

template<class T>
class Insertion : public IndelModel::Element
{
public:
	typedef T dist_type;
	typedef typename T::result_type result_type;
	
	Insertion(double dLambda = 0.0, const dist_type &dist = dist_type(), 
		const SequenceFactory &fac = SequenceFactory()) :
		m_dLambda(dLambda), rand_len(g_rng, dist), rand_seq(fac)
	{
		
	}

	virtual double Rate(const Sequence& seq) const
	{
		return m_dLambda*(seq.Length()+1);
	}

	virtual void operator()(Sequence& seq)
	{
		seq.Insert(rand_uint(seq.Length()), rand_seq(rand_len()));
	}

private:
	double m_dLambda;
	boost::variate_generator<DawgRng&, dist_type> rand_len;
	SequenceFactory rand_seq;
};

template<class T>
class Deletion : public IndelModel::Element
{
public:
	typedef T dist_type;
	typedef typename T::result_type result_type;

	Deletion(double dLambda = 0.0, const dist_type &dist = dist_type()) :
		m_dLambda(dLambda), rand_len(g_rng, dist) {	}

	virtual double Rate(const Sequence& seq) const
	{
		return m_dLambda*seq.Length();
	}
	virtual void operator()(Sequence& seq)
	{
		//do Deletion
		Sequence::pos_type uLen, uP, uB, uE;
		uLen = rand_len();
		if(uLen == Sequence::pos_type(1))
		{
			// accelerate most common type
			uB = rand_uint(seq.Length()-1);
			uE = uB+1;
		}
		else
		{
			uP = rand_uint(seq.Length()+uLen-2)+1;
			uB = std::max(uLen, uP)-uLen; 
			uE = std::min(seq.Length(), uP);
		}
		seq.Delete(uB, uE);
	}

private:
	double m_dLambda;
	boost::variate_generator<DawgRng&, dist_type> rand_len;
};

}

#endif
