#include "quaternion.H"
#include "septernion.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    quaternion q(vector(1, 0, 0), 0.7853981);
    Info<< "q " << q << endl;

    vector v(0, 1, 0);
    Info<< "v " << v << endl;

    Info<< "inv(q)*q " << inv(q)*q << endl;

    Info<< "q*quaternion(0, v)*conjugate(q) "
        << q*quaternion(0, v)*conjugate(q) << endl;

    Info<< "q.transform(v) " << q.transform(v) << endl;
    Info<< "q.R() & v " << (q.R() & v) << endl;

    Info<< "q.invTransform(v) " << q.invTransform(v) << endl;

    septernion tr(vector(0, 0.1, 0), q);
    Info<< "tr " << tr << endl;

    Info<< "inv(tr)*tr " << inv(tr)*tr << endl;

    Info<< "tr.transform(v) " << tr.transform(v) << endl;

    Info<< "(septernion(vector(0, -1, 0))*q*septernion(vector(0, 1, 0)))"
        << ".transform(v) "
        <<  (septernion(vector(0, -1, 0))
           *q
           *septernion(vector(0, 1, 0))).transform(v)
        << endl;

    return 0;
}
