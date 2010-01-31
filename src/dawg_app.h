/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 *                                                                          *
 *  This program is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 ****************************************************************************/

#ifndef DAWG_APP_H
#define DAWG_APP_H

#include <vector>
#include <utility>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/****************************************************************************
 *    class dawg_app                                                        *
 ****************************************************************************/

class dawg_app  {
public:
	dawg_app(int argc, char *argv[]);
	virtual ~dawg_app() { }
	
	virtual int run();
	
	po::options_description desc;	
	
	struct args
	{			
	// use X-Macros to specify argument variables
#	define XM(lname, sname, desc, type, def) type _V(lname) ;
#	include "dawgarg.xmh"
#	undef XM
	};
	
protected:
	args arg;
	
};

namespace boost { namespace program_options {
template<>
typed_value<bool>* value(bool* v) {
	return bool_switch(v);
}
}}

#endif
