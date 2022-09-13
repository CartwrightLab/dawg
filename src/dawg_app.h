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

#include <CLI11.hpp>
#include <utility>
#include <vector>

/****************************************************************************
 *    class dawg_app                                                        *
 ****************************************************************************/

class dawg_app {
   public:
    dawg_app(int argc, char* argv[]);
    dawg_app(const dawg_app&) = delete;             // copy constructor
    dawg_app& operator=(const dawg_app&) = delete;  // copy assignment operator
    dawg_app(dawg_app&&) = delete;                  // move constructor
    dawg_app& operator=(dawg_app&&) = delete;       // move assignment operator
    virtual ~dawg_app() = default;                  // destructor

    virtual int run();

    struct args {
        // use X-Macros to specify argument variables
#define XM(lname, sname, desc, type, def) type XV(lname);
#define XF(lname, sname, desc, type, def) type XV(lname);
#include "dawgarg.xmh"
#undef XM
#undef XF
        std::vector<std::string> input;
    };

   private:
    args arg;
    CLI::App cli_app;
    std::string runname{""};
};

#endif
