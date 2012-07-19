#pragma once
#ifndef DAWG_ZETA_H
#define DAWG_ZETA_H
/****************************************************************************
 *  Copyright (C) 2009,2012 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/math/special_functions/zeta.hpp>

namespace dawg {

inline double zeta(double z) {
	return boost::math::zeta<double>(z);
}

};

#endif

