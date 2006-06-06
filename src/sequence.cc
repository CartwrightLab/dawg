#include "sequence.h"

bool Dawg::Sequence::Replace(pos_type uPos, const residue_type &uRes)
{
	if(uPos >= vSeq.size())
		return false;
	vSeq[uPos] = uRes;
	return true;
}

bool Dawg::Sequence::Insert(pos_type uPos, const Sequence& seq)
{
	if(uPos > vSeq.size())
		return false;
	seq_type::size_type sz = seq.vSeq.size();
	vSeq.insert(vSeq.begin()+uPos, seq.vSeq.begin(), seq.vSeq.end());
	if(uPos == vSeq.size())
	{
		if(!ssAln.empty())
			ssAln.append(sz, '+');
		if(!vuSeqSections.empty())
		{
			vuSeqSections.back() += sz;
			if(!vuAlnSections.empty())
				vuAlnSections.back() += sz;
		}
	}
	else
	{
		if(!ssAln.empty())
		{
			pos_type u = SeqPosToAlnPos(uPos);
			ssAln.insert(u, sz, '+');
		}
		if(!vuSeqSections.empty())
		{
			pos_type u;
			for(u = 0; uPos > vuSeqSections[u] && u < vuSeqSections.size(); ++u)
				uPos -= vuSeqSections[u];
			vuSeqSections[u] += sz;
			if(!vuAlnSections.empty())
				vuAlnSections[u] += sz;
		}
	}
	return true;
}

bool Dawg::Sequence::Delete(pos_type uBegin, pos_type uEnd)
{
	if(uBegin > vSeq.size() || uEnd > vSeq.size() || uBegin > uEnd)
		return false;
	vSeq.erase(vSeq.begin()+uBegin, vSeq.begin()+uEnd);
	pos_type u = SeqPosToAlnPos(uBegin);
	pos_type v = SeqPosToAlnPos(uEnd);
	if(!ssAln.empty())
	{
		for(pos_type p = u; p != v; ++p)
		{
			if(ssAln[p] == '.')
				ssAln[p] = '-';
			else if(ssAln[p] == '+')
				ssAln[p] = '=';
		}
	}
	if(!vuSeqSections.empty())
	{
		for(u = 0; u < vuSeqSections.size() && uBegin > vuSeqSections[u]; ++u)
		{
			uBegin -= vuSeqSections[u];
			uEnd -= vuSeqSections[u];
		}
		if(u >= vuSeqSections.size())
			return false; // overflow
		if(uEnd <= vuSeqSections[u])
			vuSeqSections[u] -= (uEnd-uBegin);
		else
		{
			vuSeqSections[u] = uBegin;
			uEnd -= vuSeqSections[u];
			for(u += 1; u < vuSeqSections.size() && uEnd > vuSeqSections[u]; ++u)
			{
				vuSeqSections[u] = 0;
				uEnd -= vuSeqSections[u];
			}
			if(u >= vuSeqSections.size())
				return false; // overflow
			vuSeqSections[u] -= uEnd;
		}
	}
	return true;
}

//Dawg::Sequence::pos_type Dawg::Sequence::FindSection(pos_type uPos) const
//{
//	for(pos_type u = 0; u < vuSeqSections.size(); ++u)
//	{
//		if(uPos < vuSeqSections[u])
//			return u;
//		uPos -= vuSeqSections[u];
//	}
//	return static_cast<pos_type>(-1);
//}

 Dawg::Sequence::pos_type Dawg::Sequence::SeqPosToAlnPos(pos_type uPos) const
{
	for(pos_type u = 0; u < ssAln.length(); ++u)
	{
		char ch = ssAln[u];
		if((ch == '.' || ch == '+') && uPos-- == 0)
			return u;
	}
	return static_cast<pos_type>(-1);
}


bool Dawg::Sequence::CopySection(pos_type uSec, Sequence& seq) const
{
	if(uSec >= vuSeqSections.size())
		return false;
	seq.rm = rm;
	pos_type uBegin = 0;
	for(std::vector<pos_type>::const_iterator cit = vuSeqSections.begin();
		cit != vuSeqSections.begin()+uSec; ++cit)
		uBegin += *cit;
	pos_type uEnd = uBegin + vuSeqSections[uSec];
	seq.vSeq.assign(vSeq.begin()+uBegin,vSeq.begin()+uEnd);
	seq.vuSeqSections.assign(1, uEnd-uBegin);
	uBegin = 0;
	for(std::vector<pos_type>::const_iterator cit = vuAlnSections.begin();
		cit != vuAlnSections.begin()+uSec; ++cit)
		uBegin += *cit;
	uEnd = uBegin + vuAlnSections[uSec];
	seq.ssAln.assign(ssAln);
	seq.vuAlnSections.assign(1, uEnd-uBegin);
	return true;
}

bool Dawg::Sequence::AssignSection(pos_type uSec, const Sequence& seq)
{
	if(uSec >= vuSeqSections.size() || uSec >= vuAlnSections.size() )
		return false;
	pos_type uOldLen = vuSeqSections[uSec];
	pos_type uNewLen = seq.vSeq.size();
	pos_type uBegin = 0;
	for(std::vector<pos_type>::const_iterator cit = vuSeqSections.begin();
		cit != vuSeqSections.begin()+uSec; ++cit)
		uBegin += *cit;
	pos_type uEnd = uBegin + uOldLen;
	if(uOldLen == uNewLen)
	{
		for(pos_type i = 0; i < uNewLen; ++i)
			vSeq[uBegin+i] = seq.vSeq[i];
	}
	else if(uOldLen > uNewLen)
	{
		for(pos_type i = 0; i < uNewLen; ++i)
			vSeq[uBegin+i] = seq.vSeq[i];
		vSeq.erase(vSeq.begin()+uBegin+uNewLen, vSeq.begin()+uBegin+uOldLen);
	}
	else
	{
		for(pos_type i = 0; i < uOldLen; ++i)
			vSeq[uBegin+i] = seq.vSeq[i];
		vSeq.insert(vSeq.begin()+uBegin+uOldLen, seq.vSeq.begin()+uOldLen, seq.vSeq.end());
	}
	vuSeqSections[uSec] = uNewLen;

	uOldLen = vuAlnSections[uSec];
	uNewLen = seq.ssAln.size();
	uBegin = 0;
	for(std::vector<pos_type>::const_iterator cit = vuAlnSections.begin();
		cit != vuAlnSections.begin()+uSec; ++cit)
		uBegin += *cit;
	uEnd = uBegin + uOldLen;
	if(uOldLen == uNewLen)
	{
		for(pos_type i = 0; i < uNewLen; ++i)
			ssAln[uBegin+i] = seq.ssAln[i];
	}
	else if(uOldLen > uNewLen)
	{
		for(pos_type i = 0; i < uNewLen; ++i)
			ssAln[uBegin+i] = seq.ssAln[i];
		ssAln.erase(ssAln.begin()+uBegin+uNewLen, ssAln.begin()+uBegin+uOldLen);
	}
	else
	{
		for(pos_type i = 0; i < uOldLen; ++i)
			ssAln[uBegin+i] = seq.ssAln[i];
		ssAln.insert(ssAln.begin()+uBegin+uOldLen, seq.ssAln.begin()+uOldLen, seq.ssAln.end());
	}
	vuAlnSections[uSec] = uNewLen;
	return true;
}

bool Dawg::Sequence::AppendSection(const Sequence& seq)
{
	vSeq.insert(vSeq.end(), seq.vSeq.begin(), seq.vSeq.end());
	ssAln.append(seq.ssAln);
	vuSeqSections.push_back(seq.vSeq.size());
	vuAlnSections.push_back(seq.ssAln.length());
	return true;
}
