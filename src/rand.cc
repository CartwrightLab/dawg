// rand.cc - Copyright (C) 2004-2006 Reed A. Cartwright (all rights reserved)

#include "rand.h"

DawgRng g_rng;
boost::uniform_01<DawgRng, double> g_randReal01(g_rng);


