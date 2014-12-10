/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Q

Description
    Calculates and writes the second invariant of the velocity gradient tensor.

        Q = 0.5*(sqr(tr(gradU)) - tr(((gradU)&(gradU))))      [1/s^2]

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Uheader.headerOk())
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);
        volTensorField gradU(fvc::grad(U));

        volScalarField Q
        (
            IOobject
            (
                "Q",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            0.5*(sqr(tr(gradU)) - tr(((gradU)&(gradU))))
        );

        /*
        // This is a second way of calculating Q, that delivers results
        // very close, but not identical to the first approach.

        volSymmTensorField S(symm(gradU));  // symmetric part of tensor
        volTensorField W(skew(gradU));      // anti-symmetric part

        volScalarField SS(S && S);
        volScalarField WW(W && W);

        volScalarField Q
        (
            IOobject
            (
                "Q",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            0.5*(WW - SS)
        );
        */

        Info<< "mag(Q) max/min : "
            << max(Q).value() << " "
            << min(Q).value() << endl;

        if (writeResults)
        {
            Q.write();
        }
    }
    else
    {
        Info<< "    No U" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
