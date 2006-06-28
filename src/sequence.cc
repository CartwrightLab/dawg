//#pragma warning(disable: 4127 4512)

#include "sequence.h"

void Dawg::Sequence::Clear()
{
	m_vBases.clear();
	m_vRates.clear();
	m_vSeqSections.clear();
	m_vAln.clear();
	m_vAlnSections.clear();
}

bool Dawg::Sequence::Replace(pos_type uPos, const base_type &base, const rate_type &rate)
{
	if(uPos >= m_vBases.size())
		return false;
	m_vBases[uPos] = base;
	if(uPos < m_vRates.size())
		m_vRates[uPos] = rate;
	return true;
}

bool Dawg::Sequence::Insert(pos_type uPos, const Sequence& seq)
{
	if(uPos > m_vBases.size())
		return false;
	pos_type sz = seq.m_vBases.size();
	if(!m_vAln.empty())
	{
		if(uPos == Length())
			m_vAln.insert(m_vAln.end(), sz, '+');
		else
		{
			pos_type u = SeqPosToAlnPos(uPos);
			m_vAln.insert(m_vAln.begin()+u, sz, '+');			
		}
	}
	if(!m_vSeqSections.empty())
	{
		if(uPos == Length())
		{
			m_vSeqSections.back() += sz;
			if(!m_vAlnSections.empty())
				m_vAlnSections.back() += sz;
		}
		else
		{
			pos_type u;
			for(u = 0; uPos > m_vSeqSections[u] && u < m_vSeqSections.size(); ++u)
				uPos -= m_vSeqSections[u];
			m_vSeqSections[u] += sz;
			if(!m_vAlnSections.empty())
				m_vAlnSections[u] += sz;
		}
	}
	m_vBases.insert(m_vBases.begin()+uPos, seq.m_vBases.begin(), seq.m_vBases.end());
	if(!m_vRates.empty())
		m_vRates.insert(m_vRates.begin()+uPos, seq.m_vRates.begin(), seq.m_vRates.end());
	return true;
}

bool Dawg::Sequence::Delete(pos_type uBegin, pos_type uEnd)
{
	if(uBegin > m_vBases.size() || uEnd > m_vBases.size() || uBegin > uEnd)
		return false;
	m_vBases.erase(m_vBases.begin()+uBegin, m_vBases.begin()+uEnd);
	m_vRates.erase(m_vRates.begin()+uBegin, m_vRates.begin()+uEnd);

	pos_type u = SeqPosToAlnPos(uBegin);
	pos_type v = SeqPosToAlnPos(uEnd);
	if(!m_vAln.empty())
	{
		for(pos_type p = u; p != v; ++p)
		{
			if(m_vAln[p] == '.')
				m_vAln[p] = '-';
			else if(m_vAln[p] == '+')
				m_vAln[p] = '=';
		}
	}
	if(!m_vSeqSections.empty())
	{
		for(u = 0; u < m_vSeqSections.size() && uBegin > m_vSeqSections[u]; ++u)
		{
			uBegin -= m_vSeqSections[u];
			uEnd -= m_vSeqSections[u];
		}
		if(u >= m_vSeqSections.size())
			return false; // overflow
		if(uEnd <= m_vSeqSections[u])
			m_vSeqSections[u] -= (uEnd-uBegin);
		else
		{
			m_vSeqSections[u] = uBegin;
			uEnd -= m_vSeqSections[u];
			for(u += 1; u < m_vSeqSections.size() && uEnd > m_vSeqSections[u]; ++u)
			{
				m_vSeqSections[u] = 0;
				uEnd -= m_vSeqSections[u];
			}
			if(u >= m_vSeqSections.size())
				return false; // overflow
			m_vSeqSections[u] -= uEnd;
		}
	}
	return true;
}

Dawg::Sequence::pos_type Dawg::Sequence::SeqPosToAlnPos(pos_type uPos) const
{
	for(pos_type u = 0; u < m_vAln.size(); ++u)
	{
		char ch = m_vAln[u];
		if((ch == '.' || ch == '+') && uPos-- == 0)
			return u;
	}
	return static_cast<pos_type>(-1);
}

bool Dawg::Sequence::AppendSection(const Sequence& seq)
{
	if(m_vSeqSections.empty())
		m_vSeqSections.push_back(m_vBases.size());
	if(!m_vAln.empty())
	{
		m_vAln.insert(m_vAln.end(), seq.m_vAln.begin(), seq.m_vAln.end());
		if(m_vAlnSections.empty())
			m_vAlnSections.push_back(m_vAln.size());
	}
	m_vBases.insert(m_vBases.end(), seq.m_vBases.begin(), seq.m_vBases.end());
	if(!m_vRates.empty())
		m_vRates.insert(m_vRates.end(), seq.m_vRates.begin(), seq.m_vRates.end());
	m_vSeqSections.push_back(seq.m_vBases.size());
	if(!m_vAlnSections.empty())
		m_vAlnSections.push_back(seq.m_vAln.size());
	return true;
}

bool Dawg::Sequence::ReadSection(pos_type uSec, Sequence &seq)
{
	if(uSec >= m_vSeqSections.size())
		return false;
	pos_type uBegin = 0;
	for(pos_type u = 0; u < uSec; ++u)
		uBegin += m_vSeqSections[u];
	pos_type uEnd = uBegin + m_vSeqSections[uSec];
	seq.m_vBases.assign(m_vBases.begin()+uBegin, m_vBases.begin()+uEnd);
	if(m_vRates.empty())
		seq.m_vRates.clear();
	else
		seq.m_vRates.assign(m_vRates.begin()+uBegin, m_vRates.begin()+uEnd);
	seq.m_vSeqSections.assign(1, uEnd-uBegin);
	if(m_vAln.empty())
	{
		seq.m_vAln.clear();
		seq.m_vAlnSections.clear();
	}
	else
	{
		uBegin = 0;
		for(pos_type u = 0; u < uSec; ++u)
			uBegin += m_vAlnSections[u];
		uEnd = uBegin + m_vAlnSections[uSec];
		seq.m_vAln.assign(m_vAln.begin()+uBegin, m_vAln.begin()+uEnd);
		seq.m_vAlnSections.assign(1, uEnd-uBegin);
	}

}

void Dawg::SequenceFactory::operator()(Sequence& seq, Sequence::pos_type uLen)
{
	seq.m_vBases.resize(uLen);
	std::generate(seq.m_vBases.begin(), seq.m_vBases.end(), rand_base);
	seq.m_vRates.resize(uLen);
	std::generate(seq.m_vRates.begin(), seq.m_vRates.end(), rand_rate);
	seq.m_vAln.assign(uLen, '.');
	seq.m_vAlnSections.assign(1, uLen);
	seq.m_vSeqSections.assign(1, uLen);
}

Dawg::Sequence Dawg::SequenceFactory::operator ()(Sequence::pos_type uLen)
{
	Sequence seq;
	(*this)(seq, uLen);
	return seq;
}



