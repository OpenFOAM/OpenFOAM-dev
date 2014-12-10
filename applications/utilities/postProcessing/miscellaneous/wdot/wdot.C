/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    wdot

Description
    Calculates and writes wdot for each time.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        mesh.readUpdate();

        volScalarField mgb
        (
            IOobject
            (
                "mgb",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        volScalarField Xi
        (
            IOobject
            (
                "Xi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        volScalarField St
        (
            IOobject
            (
                "St",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            Xi*Su
        );

        St.write();

        volScalarField wdot
        (
            IOobject
            (
                "wdot",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
           St*mgb
        );

        wdot.write();
    }

    return 0;
}


// ************************************************************************* //
