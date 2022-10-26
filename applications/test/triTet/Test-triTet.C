#include "point.H"
#include "triangle.H"
#include "tetrahedron.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    const triangle<point, point> tri
    (
        vector(0, 0, 0),
        vector(1, 0, 0),
        vector(1, 1, 0)
    );

    const Tuple2<point, scalar> triCircle = tri.circumCircle();

    Info<< "tri circumCentre = " << triCircle.first() << endl
        << "tri circumRadius = " << triCircle.second() << endl
        << "     tri quality = " << tri.quality() << endl;

    const tetrahedron<point, point> tet
    (
        vector(1, 0, 0),
        vector(0, 1, 0),
        vector(0, 0, 1),
        vector(0.5773502, 0.5773502, 0.5773502)
    );

    const Tuple2<point, scalar> tetSphere = tet.circumSphere();

    Info<< "tet circumCentre = " << tetSphere.first() << endl
        << "tet circumRadius = " << tetSphere.second() << endl
        << "     tet quality = " << tet.quality() << endl;

    return 0;
}
