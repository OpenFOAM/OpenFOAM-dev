/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Test-tetTetOverlap

Description
    Overlap volume of two tets

\*---------------------------------------------------------------------------*/

#include "tetrahedron.H"
#include "OFstream.H"
#include "meshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeOBJ
(
    Ostream& os,
    label& vertI,
    const tetPoints& tet
)
{
    forAll(tet, fp)
    {
        meshTools::writeOBJ(os, tet[fp]);
    }
    os  << "l " << vertI+1 << ' ' << vertI+2 << nl
        << "l " << vertI+1 << ' ' << vertI+3 << nl
        << "l " << vertI+1 << ' ' << vertI+4 << nl
        << "l " << vertI+2 << ' ' << vertI+3 << nl
        << "l " << vertI+2 << ' ' << vertI+4 << nl
        << "l " << vertI+3 << ' ' << vertI+4 << nl;
    vertI += 4;
}


int main(int argc, char *argv[])
{
    tetPoints A
    (
        point(0, 0, 0),
        point(1, 0, 0),
        point(1, 1, 0),
        point(1, 1, 1)
    );
    const tetPointRef tetA = A.tet();

    tetPoints B
    (
        point(0.1, 0.1, 0.1),
        point(1.1, 0.1, 0.1),
        point(1.1, 1.1, 0.1),
        point(1.1, 1.1, 1.1)
    );
    const tetPointRef tetB = B.tet();


    tetPointRef::tetIntersectionList insideTets;
    label nInside = 0;
    tetPointRef::tetIntersectionList outsideTets;
    label nOutside = 0;

    tetA.tetOverlap
    (
        tetB,
        insideTets,
        nInside,
        outsideTets,
        nOutside
    );


    // Dump to file
    // ~~~~~~~~~~~~

    {
        OFstream str("tetA.obj");
        Info<< "Writing A to " << str.name() << endl;
        label vertI = 0;
        writeOBJ(str, vertI, A);
    }
    {
        OFstream str("tetB.obj");
        Info<< "Writing B to " << str.name() << endl;
        label vertI = 0;
        writeOBJ(str, vertI, B);
    }
    {
        OFstream str("inside.obj");
        Info<< "Writing parts of A inside B to " << str.name() << endl;
        label vertI = 0;
        for (label i = 0; i < nInside; ++i)
        {
            writeOBJ(str, vertI, insideTets[i]);
        }
    }
    {
        OFstream str("outside.obj");
        Info<< "Writing parts of A outside B to " << str.name() << endl;
        label vertI = 0;
        for (label i = 0; i < nOutside; ++i)
        {
            writeOBJ(str, vertI, outsideTets[i]);
        }
    }


    // Check
    // ~~~~~

    Info<< "Vol A:" << tetA.mag() << endl;

    scalar volInside = 0;
    for (label i = 0; i < nInside; ++i)
    {
        volInside += insideTets[i].tet().mag();
    }
    Info<< "Vol A inside B:" << volInside << endl;

    scalar volOutside = 0;
    for (label i = 0; i < nOutside; ++i)
    {
        volOutside += outsideTets[i].tet().mag();
    }
    Info<< "Vol A outside B:" << volOutside << endl;

    Info<< "Sum inside and outside:" << volInside+volOutside << endl;

    if (mag(volInside+volOutside-tetA.mag()) > SMALL)
    {
        FatalErrorIn("Test-tetetOverlap")
            << "Tet volumes do not sum up to input tet."
            << exit(FatalError);
    }

    return 0;
}


// ************************************************************************* //
