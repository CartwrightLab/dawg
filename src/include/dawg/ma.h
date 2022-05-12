#pragma once
#ifndef DAWG_MA_H
#define DAWG_MA_H
/****************************************************************************
 *  Copyright (C) 2009 Reed A. Cartwright, PhD <reed@scit.us>               *
 ****************************************************************************/

#include <boost/algorithm/string/case_conv.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "dawg/trick.h"
#include "dawg/utils/vecio.h"

namespace dawg {

// dawg::ma is a "model argument" structure
struct ma {
#define XM(name, type, def, desc) type XV(name);
#include "dawg/details/dawgma.xmh"
#undef XM
    std::string name;

    ma(const std::string &_n = std::string())
        :
#define XM(name, type, def, desc) XV(name)(def),
#include "dawg/details/dawgma.xmh"
#undef XM
          name(_n) {
    }

    static bool from_trick(const dawg::trick &trk, std::vector<dawg::ma> &v);

    void read_section(const trick::data_type::value_type &sec);

    template <class CharType, class CharTrait>
    static void help(std::basic_ostream<CharType, CharTrait> &o);

   private:
};

template <class CharType, class CharTrait>
inline std::basic_ostream<CharType, CharTrait> &operator<<(
    std::basic_ostream<CharType, CharTrait> &o, const ma &a) {
    if(!o.good()) return o;

    o << set_open('\x7f') << set_close('\x7f') << set_delimiter(',');

    o << "[[ " << a.name << " ]]" << std::endl;
#define XM(name, type, def, desc) \
    o << XP(name) " = " << a.XV(name) << std::endl;
#include "dawg/details/dawgma.xmh"
#undef XM

    return o;
}

namespace details {
inline std::string ma_help_name(const char *cs) {
    std::string ret(cs);
    typedef boost::iterator_range<std::string::iterator> range;
    range r(ret.begin(), ret.begin() + 1);
    boost::to_upper(r);
    std::size_t pos = 0;
    while((pos = ret.find('.', pos + 1)) != std::string::npos) {
        r = range(ret.begin() + pos + 1, ret.begin() + pos + 2);
        boost::to_upper(r);
    }
    return ret;
}
};  // namespace details

template <class CharType, class CharTrait>
void ma::help(std::basic_ostream<CharType, CharTrait> &o) {
    o << "[REGULAR PARAMETERS]" << std::endl;
#define XM(name, type, def, desc) \
    o << details::ma_help_name(XP(name)) << " - " << (desc) << std::endl;
#include "dawg/details/dawgma.xmh"
#undef XM
    o << "\n[GLOBAL PARAMETERS]" << std::endl;
#define XM(name, type, def, desc) \
    o << details::ma_help_name(XP(name)) << " - " << (desc) << std::endl;
#include "dawg/details/global.xmh"
#undef XM
}

}  // namespace dawg
#endif  // DAWG_MA_H
