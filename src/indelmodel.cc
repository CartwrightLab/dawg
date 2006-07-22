#pragma warning(disable: 4512)

#include <algorithm>
#include "indelmodel.h"

#include "sequence.h"

void Dawg::IndelModel::operator()(Sequence& seq, double dTime)
{
	if(m_elements.size() == 0 || dTime <= 0.0)
		return;
	std::vector<double> dRates(m_elements.size());
	double dSum, dW;
	while(1)
	{
		dSum = 0.0;
		for(std::vector<double>::size_type u = 0; u < dRates.size(); ++u )
			dRates[u] = (dSum += m_elements[u].Rate(seq));
		dTime -= Waiting(dSum);
		if(dTime <= 0.0)
			break;
		dW = Which(dSum);
		for(std::vector<double>::size_type u = 0; u < dRates.size(); ++u )
		{
			if( dW < dRates[u])
			{
				m_elements[u](seq);
				break;
			}
		}
	}		
}

