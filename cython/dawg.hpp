#ifndef DAWG_HPP
#define DAWG_HPP

#include <string>
#include <list>

namespace dawg {

class Dawg
{
public:
    explicit Dawg();
    explicit Dawg(const std::list<std::string>& args);
    void run();
private:

}; // class Dawg

} // namespace dawg

#endif
