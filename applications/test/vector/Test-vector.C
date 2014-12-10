#include "vector.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    Info<< vector::zero << endl
        << vector::one << endl
        << vector::dim << endl
        << vector::rank << endl;

    vector d(0.5, 0.5, 0.5);
    d /= mag(d);

    vector dSmall = (1e-100)*d;
    dSmall /= mag(dSmall);

    Info<< (dSmall - d) << endl;

    d *= 4.0;

    Info<< d << endl;

    Info<< d + d << endl;

    Info<< magSqr(d) << endl;

    vector d2(0.5, 0.51, -0.5);
    Info<< cmptMax(d2) << " "
        << cmptSum(d2) << " "
        << cmptProduct(d2) << " "
        << cmptMag(d2)
        << endl;
    Info<< min(d, d2) << endl;
    return 0;
}
