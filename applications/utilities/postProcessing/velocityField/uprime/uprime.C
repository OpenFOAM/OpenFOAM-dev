/*---------------------------------------------------------------------------* \
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
    uprime

Description
    Calculates and writes the scalar field of uprime (sqrt(2/3 k)).

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject kheader
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (kheader.headerOk())
    {
        Info<< "    Reading k" << endl;
        volScalarField k(kheader, mesh);

        Info<< "    Calculating uprime" << endl;
        volScalarField uprime
        (
            IOobject
            (
                "uprime",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            sqrt((2.0/3.0)*k)
        );

        Info<< "uprime max/min : "
            << max(uprime).value() << " "
            << min(uprime).value() << endl;

        if (writeResults)
        {
            uprime.write();
        }
    }
    else
    {
        Info<< "    No k" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
