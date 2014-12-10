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
    pPrime2

Description
    Calculates and writes the scalar field of pPrime2 (sqr(p - pMean)) at
    each time

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

    runTime.setTime(timeDirs.last(), timeDirs.size()-1);

    volScalarField pMean
    (
        IOobject
        (
            "pMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    );

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject pheader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check p exists
        if (pheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading p" << endl;
            volScalarField p(pheader, mesh);

            Info<< "    Calculating pPrime2" << endl;
            volScalarField pPrime2
            (
                IOobject
                (
                    "pPrime2",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                sqr(p - pMean)
            );
            pPrime2.write();
        }
        else
        {
            Info<< "    No p" << endl;
        }

        Info<< endl;
    }

    return 0;
}


// ************************************************************************* //
