#ifndef DAWG_SEQUENCE_H
#define DAWG_SEQUENCE_H

#include <boost/numeric/ublas/matrix.hpp>
#include "rand.h"

namespace ublas = boost::numeric::ublas;

namespace Dawg {

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
};

class SequenceFactory
{
public:
	SequenceFactory() : rand_rate(g_rng, gammaiota_distribution<>()) ,
		rand_base(g_rng, discrete_distribution<Sequence::pos_type>()) { }

	void operator()(Sequence& seq, Sequence::pos_type uLen);

protected:
	boost::variate_generator<DawgRng&, gammaiota_distribution<> > rand_rate;
	boost::variate_generator<DawgRng&, discrete_distribution<Sequence::pos_type> > rand_base;
};

}; // namespace Dawg

#endif
