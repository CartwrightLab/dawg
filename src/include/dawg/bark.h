#pragma once
#ifndef DAWG_BARK_H
#define DAWG_BARK_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/exception/all.hpp>

namespace dawg {

struct bark: virtual boost::exception, virtual std::exception { };

} // namespace dawg


#endif

