#include "dawg.hpp"

#include "../src/dawg_app.h"
#include <dawg/ma.h>
#include <dawg/trick.h>
#include <dawg/global.h>
#include <dawg/output.h>

#include <iostream>

dawg::Dawg::Dawg()
{
    using namespace std;
    cout << "Hello PyDawg: " << __FILE__ << ", " << __LINE__ << endl;
}

dawg::Dawg::Dawg(const std::string& s)
{
    using namespace std;
    cout << "Hello PyDawg: " << __FILE__ << ", " << __LINE__ << endl;
    cout << "s: " << s << "\n";
}

void dawg::Dawg::run()
{
    using namespace std;
    cout << "Hello PyDawg" << __FILE__ << ", " << __LINE__ << endl;

}
