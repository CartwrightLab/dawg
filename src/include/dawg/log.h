#pragma once
#ifndef DAWG_LOG_H
#define DAWG_LOG_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <iostream>
#include <iomanip>

#ifdef NDEBUG
#	define DAWG_ERROR(err_msg) ((std::cerr << "ERROR: " << err_msg << std::endl), false)
#	define DAWG_WARN(warn_msg) ((std::cerr << "WARNING: " << warn_msg << std::endl), false)
#else
#	define DAWG_ERROR(err_msg) ((std::cerr << "ERROR: " << err_msg << std::endl \
	<< "(" << __FILE__ << ":" << __LINE__ << ")" << std::endl ), false)
#	define DAWG_WARN(warn_msg) ((std::cerr << "WARNING: " << warn_msg << std::endl \
	<< "(" << __FILE__ << ":" << __LINE__ << ")" << std::endl ), false)
#endif

#endif /* DAWG_LOG_H */
