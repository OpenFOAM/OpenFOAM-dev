#include "baseClass.H"

using namespace Foam;

int main(int argc, char* argv[]) {

    auto ptr = Foam::baseClass<int, float>::New("childClass", 5.0);
    ptr->correct();
    auto ptr2 = Foam::baseClass<float, double>::New("childClass", 8.0);
    ptr2->correct();
    return 0;
}
