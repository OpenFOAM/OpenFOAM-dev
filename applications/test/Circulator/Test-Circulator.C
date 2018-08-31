/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-circulator

Description

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "ListOps.H"
#include "face.H"
#include "Circulator.H"
#include "ConstCirculator.H"


using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Test the implementation of a circular iterator" << nl << endl;

    Info<< "Test const circulator. First go forwards, then backwards."
        << nl << endl;

    face f(identity(4));

    ConstCirculator<face> cStart(f);

    if (cStart.size()) do
    {
        Info<< "Iterate forwards over face (prev/curr/next) : "
            << cStart.prev() << " / " << cStart() << " / " << cStart.next()
            << endl;

    } while (cStart.circulate(CirculatorBase::direction::clockwise));

    if (cStart.size()) do
    {
        Info<< "Iterate backwards over face : " << cStart() << endl;

    } while (cStart.circulate(CirculatorBase::direction::anticlockwise));


    Info<< nl << nl << "Test non-const circulator" << nl << endl;

    Circulator<face> cStart2(f);

    Info<< "Face before : " << f << endl;

    if (cStart2.size()) do
    {
        Info<< "Iterate forwards over face (prev/curr/next) : "
            << cStart2.prev() << " / " << cStart2() << " / " << cStart2.next()
            << endl;

    } while (cStart2.circulate(CirculatorBase::direction::clockwise));

    if (cStart2.size()) do
    {
        Info<< "Iterate forwards over face, adding 1 to each element : "
            << cStart2();

        cStart2() += 1;

        Info<< " -> " << cStart2() << endl;
    } while (cStart2.circulate(CirculatorBase::direction::clockwise));

    Info<< "Face after : " << f << endl;


    Info<< nl << nl << "Compare two faces: " << endl;
    face a(identity(5));
    Info<< "Compare " << a << " and " << a << " Match = " << face::compare(a, a)
        << endl;

    face b(reverseList(a));
    Info<< "Compare " << a << " and " << b << " Match = " << face::compare(a, b)
        << endl;

    face c(a);
    c[4] = 3;
    Info<< "Compare " << a << " and " << c << " Match = " << face::compare(a, c)
        << endl;

    face d(rotateList(a, 2));
    Info<< "Compare " << a << " and " << d << " Match = " << face::compare(a, d)
        << endl;

    face g(labelList(5, 1));
    face h(g);
    Info<< "Compare " << g << " and " << h << " Match = " << face::compare(g, h)
        << endl;

    g[0] = 2;
    h[3] = 2;
    Info<< "Compare " << g << " and " << h << " Match = " << face::compare(g, h)
        << endl;

    g[4] = 3;
    h[4] = 3;
    Info<< "Compare " << g << " and " << h << " Match = " << face::compare(g, h)
        << endl;

    face face1(identity(1));
    Info<< "Compare " << face1 << " and " << face1
        << " Match = " << face::compare(face1, face1) << endl;

    face face2(identity(1)+1);
    Info<< "Compare " << face1 << " and " << face2
        << " Match = " << face::compare(face1, face2) << endl;

    Info<< nl << nl << "Zero face" << nl << endl;

    face fZero;
    Circulator<face> cZero(fZero);

    if (cZero.size()) do
    {
        Info<< "Iterate forwards over face : " << cZero() << endl;

    } while (cZero.circulate(CirculatorBase::direction::clockwise));

    fZero = face(identity(5));

    // circulator was invalidated so reset
    cZero = Circulator<face>(fZero);

    do
    {
        Info<< "Iterate forwards over face : " << cZero() << endl;

    } while (cZero.circulate(CirculatorBase::direction::clockwise));


    Info<< nl << nl << "Simultaneously go forwards/backwards over face " << f
        << nl << endl;

    ConstCirculator<face> circForward(f);
    ConstCirculator<face> circBackward(f);

    if (circForward.size() && circBackward.size()) do
    {
        Info<< "Iterate over face forwards : " << circForward()
            << ", backwards : " << circBackward() << endl;
    }
    while
    (
        circForward.circulate(CirculatorBase::direction::clockwise),
        circBackward.circulate(CirculatorBase::direction::anticlockwise)
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
