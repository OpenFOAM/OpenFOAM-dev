#include "tensor2D.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    vector2D v1(1, 2), v2(3, 4);
    tensor2D t = v1*v2;

    Info<< "v1(1, 2)*v2(3, 4) = " << t << endl;

    return 0;
}
