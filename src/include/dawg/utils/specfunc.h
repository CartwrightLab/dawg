#pragma once
#ifndef DAWG_SPECFUNC_H
#define DAWG_SPECFUNC_H
/****************************************************************************
 *  Copyright (C) 2009,2012,2013 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/beta.hpp>

namespace dawg {

inline double zeta(double z) {
	return boost::math::zeta<double>(z);
}

inline double beta(double a, double b) {
	return boost::math::beta<double,double>(a,b);
}

};

#endif

