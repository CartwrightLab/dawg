#ifndef DAWG_SEQUENCE_H
#define DAWG_SEQUENCE_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "rand.h"

namespace ublas = boost::numeric::ublas;

namespace Dawg {

class SequenceFactory;

class Sequence
{
public:
	typedef ublas::matrix<double>::size_type base_type;
	typedef std::pair<base_type, double> residue_type;
	typedef std::vector<residue_type> seq_type;
	typedef seq_type::size_type pos_type;
	enum RateMode {rmHomo, rmHetro};

	pos_type Length() const { return vSeq.size(); }
	const residue_type& Residue(pos_type uPos) const { return vSeq.at(uPos); }
	residue_type& Residue(pos_type uPos) { return vSeq.at(uPos); }
	base_type Base(pos_type uPos) const { return Residue(uPos).first; }
	double Rate(pos_type uPos) const { return (rm == rmHetro) ? Residue(uPos).second : 1.0; }

	bool Replace(pos_type uPos, const residue_type &uRes);
	bool Insert(pos_type uPos, const Sequence& seq);
	bool Delete(pos_type uBegin, pos_type uEnd);

	void Clear();

	bool CopySection(pos_type uSec, Sequence& seq) const;
	bool AssignSection(pos_type uSec, const Sequence& seq);
	bool AppendSection(const Sequence& seq);
	std::vector<pos_type>::size_type SectionCount() const { return vuSeqSections.size(); }

	pos_type FindSection(pos_type uPos) const;
	pos_type SeqPosToAlnPos(pos_type uPos) const;

protected:
	RateMode rm;
	seq_type vSeq;
	std::string ssAln;
	std::vector<pos_type> vuSeqSections;
	std::vector<pos_type> vuAlnSections;

	friend SequenceFactory;
};

class SequenceFactory
{
public:
	typedef discrete_distribution<Sequence::pos_type> base_dist;
	typedef gammaiota_distribution<> rate_dist;

	SequenceFactory(const base_dist &bd = base_dist(), const rate_dist &rd = rate_dist()) : 
		rand_base(g_rng, bd), rand_rate(g_rng, rd)
		 { }

	void operator()(Sequence& seq, Sequence::pos_type uLen);
	Sequence operator()(Sequence::pos_type uLen);

	const SequenceFactory& operator=(const SequenceFactory& right)
	{
		if(&right == this)
			return *this;
		rand_base.distribution() = right.rand_base.distribution();
		rand_rate.distribution() = right.rand_rate.distribution();
		*this;
	}

protected:
	boost::variate_generator<DawgRng&, base_dist > rand_base;
	boost::variate_generator<DawgRng&, rate_dist > rand_rate;
};

class GillespieProcessor
{
public:
	class Element
	{
	public:
		virtual double Rate(const Sequence& /*seq*/) const { return 0.0;}
		virtual void operator()(Sequence& /*seq*/) { }
		void Process(Sequence & seq) { (*this)(seq); }
	};
	
	void AddElement(Element* p) { m_elements.push_back(p); }
	
	void operator()(Sequence& seq, double dTime);

private:
	double Waiting(double d) { return rand_exp(d); }
	double Which(double d) { return rand_real(0.0, d); }

	boost::ptr_vector<Element> m_elements;
};

}; // namespace Dawg

#endif
