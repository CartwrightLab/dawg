#ifndef DAWG_SEQUENCE_H
#define DAWG_SEQUENCE_H

#include <boost/numeric/ublas/matrix.hpp>
#include "rand.h"

namespace ublas = boost::numeric::ublas;

namespace Dawg {

class SequenceFactory;

class Sequence
{
public:
	typedef ublas::matrix<double>::size_type base_type;
	typedef double rate_type;
	typedef std::vector<base_type> base_seq_type;
	typedef std::vector<rate_type> rate_seq_type;
	typedef std::vector<char> aln_seq_type;
	typedef base_seq_type::size_type pos_type;
	typedef std::vector<pos_type> section_seq_type;

	pos_type Length() const { return m_vBases.size(); }
	base_type Base(pos_type uPos) const { return m_vBases[uPos]; }
	rate_type Rate(pos_type uPos) const { return (m_vRates.size() ? m_vRates[uPos] : 1.0); }
	
	bool Replace(pos_type uPos, const base_type &base, const rate_type &rate);
	bool Insert(pos_type uPos, const Sequence& seq);
	bool Delete(pos_type uBegin, pos_type uEnd);

	const aln_seq_type& Aln() const { return m_vAln; }

	void Clear();

	// Add seq as a new section
	bool AppendSection(const Sequence& seq);
	// Extract section as sequence
	bool ReadSection(pos_type uSec, Sequence & seq);
	
	section_seq_type::size_type SectionCount() const { return m_vSeqSections.size(); }

	pos_type FindSection(pos_type uPos) const;
	pos_type SeqPosToAlnPos(pos_type uPos) const;

protected:
	base_seq_type m_vBases;
	rate_seq_type m_vRates;
	section_seq_type m_vSeqSections;

	aln_seq_type  m_vAln;
	section_seq_type m_vAlnSections;

	friend SequenceFactory;
};

class SequenceFactory
{
public:
	enum Flags {
		FlagDefault    = 0,
		FlagRateHetero = 1,
		FlagAlignment  = 2,
		FlagSections   = 4
	};

	typedef discrete_distribution<Sequence::pos_type> base_dist;
	typedef gammaiota_distribution<> rate_dist;

	SequenceFactory(const base_dist &bd = base_dist(), const rate_dist &rd = rate_dist()) : 
		rand_base(g_rng, bd), rand_rate(g_rng, rd)
		 { }
	
	bool Create(const std::vector<double> &vdFreqs, double dGamma, double dIota)
	{
		rand_base.distribution() = base_dist(vdFreqs.begin(), vdFreqs.end());
		rand_rate.distribution() = rate_dist(dGamma, dIota);
		return true;
	}


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

}; // namespace Dawg

#endif
