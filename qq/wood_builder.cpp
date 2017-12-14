#include "wood_parse.h"
#include "wood.h"

#include <string>

int main(void) {
    using namespace std;

    static string tree1 = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);";
    static string tree2 = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);";
    static string tree3 = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);";
    static string tree4 = "((Man:0.1,Monkey:0.1):0.2,Dawg:0.25);";

    dawg::wood usertree;
    if (!usertree.parse(tree1.begin(), tree1.end())) {
        cout << "parse failed for tree: " << tree1 << endl;
    }

    auto tree1Parsed = usertree.data();
    for (auto &&node : tree1Parsed) {
        cout << "node: " << node << endl;
    }

    return 0;
}
