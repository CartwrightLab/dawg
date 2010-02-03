#pragma once
#ifndef DAWG_ZETA_H
#define DAWG_ZETA_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <gsl/gsl_sf_zeta.h>

namespace dawg {

double zeta(double z) {
	return gsl_sf_zeta(z);
}


};

#endif

