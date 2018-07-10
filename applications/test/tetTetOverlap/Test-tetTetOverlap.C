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
    Test-tetTetOverlap

Description
    Overlap volume of two tets

\*---------------------------------------------------------------------------*/

#include "tetPointRef.H"
#include "OFstream.H"
#include "meshTools.H"
#include "cut.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeOBJ
(
    Ostream& os,
    label& vertI,
    const FixedList<point, 4>& tet
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


tetPointRef makeTetPointRef(const FixedList<point, 4>& p)
{
    return tetPointRef(p[0], p[1], p[2], p[3]);
}


int main(int argc, char *argv[])
{
    // Tets to test
    FixedList<point, 4> tetA
    ({
        point(0, 0, 0),
        point(1, 0, 0),
        point(1, 1, 0),
        point(1, 1, 1)
    });
    FixedList<point, 4> tetB
    ({
        point(0.1, 0.1, 0.1),
        point(1.1, 0.1, 0.1),
        point(1.1, 1.1, 0.1),
        point(1.1, 1.1, 1.1)
    });


    // Do intersection
    typedef DynamicList<FixedList<point, 4>> tetList;
    tetList tetsIn1, tetsIn2, tetsOut;
    cut::appendOp<tetList> tetOpIn1(tetsIn1);
    cut::appendOp<tetList> tetOpIn2(tetsIn2);
    cut::appendOp<tetList> tetOpOut(tetsOut);

    const plane p0(tetB[1], tetB[3], tetB[2]);
    tetsIn1.clear();
    tetCut(tetA, p0, tetOpIn1, tetOpOut);

    const plane p1(tetB[0], tetB[2], tetB[3]);
    tetsIn2.clear();
    forAll(tetsIn1, i)
    {
        tetCut(tetsIn1[i], p1, tetOpIn2, tetOpOut);
    }

    const plane p2(tetB[0], tetB[3], tetB[1]);
    tetsIn1.clear();
    forAll(tetsIn2, i)
    {
        tetCut(tetsIn2[i], p2, tetOpIn1, tetOpOut);
    }

    const plane p3(tetB[0], tetB[1], tetB[2]);
    tetsIn2.clear();
    forAll(tetsIn1, i)
    {
        tetCut(tetsIn1[i], p3, tetOpIn2, tetOpOut);
    }

    const tetList& tetsIn = tetsIn2;


    // Dump to file
    {
        OFstream str("A.obj");
        Info<< "Writing A to " << str.name() << endl;
        label vertI = 0;
        writeOBJ(str, vertI, tetA);
    }
    {
        OFstream str("B.obj");
        Info<< "Writing B to " << str.name() << endl;
        label vertI = 0;
        writeOBJ(str, vertI, tetB);
    }
    {
        OFstream str("AInB.obj");
        Info<< "Writing parts of B inside A to " << str.name() << endl;
        label vertI = 0;
        forAll(tetsIn, i)
        {
            writeOBJ(str, vertI, tetsIn[i]);
        }
    }
    {
        OFstream str("AOutB.obj");
        Info<< "Writing parts of B inside A to " << str.name() << endl;
        label vertI = 0;
        forAll(tetsOut, i)
        {
            writeOBJ(str, vertI, tetsOut[i]);
        }
    }


    // Check the volumes
    Info<< "Vol A: " << makeTetPointRef(tetA).mag() << endl;

    scalar volIn = 0;
    forAll(tetsIn, i)
    {
        volIn += makeTetPointRef(tetsIn[i]).mag();
    }
    Info<< "Vol A inside B: " << volIn << endl;

    scalar volOut = 0;
    forAll(tetsOut, i)
    {
        volOut += makeTetPointRef(tetsOut[i]).mag();
    }
    Info<< "Vol A outside B: " << volOut << endl;

    Info<< "Sum inside and outside: " << volIn + volOut << endl;

    if (mag(volIn + volOut - makeTetPointRef(tetA).mag()) > small)
    {
        FatalErrorInFunction
            << "Tet volumes do not sum up to input tet."
            << exit(FatalError);
    }


    return 0;
}


// ************************************************************************* //
