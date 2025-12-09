#include "../src/MLFMA.h"
#include "../src/math.h"
#include "../src/node.h"

using namespace std;

void testINodeFuncs() {

    cout << "\nTesting INodeDistances()...\n";
    
    auto dists = Math::getINodeDistances();

    for (int i = 0; i < dists.size(); ++i)
        cout << i << ' ' << dists[i] << '\n';

    cout << "\nTesting INodeDirections()...\n";

    auto dirs = Math::getINodeDirections();

    for (int i = 0; i < dirs.size(); ++i)
        cout << i << ' ' << dirs[i] << '\n';

    // auto dr = dirs[254]-dirs[310];
    // cout << dr << ' ' << dr.norm() << ' ' << (dr.norm() < 1.0E-2) << '\n';
}

int main() {

    testINodeFuncs();

    return 0;
}