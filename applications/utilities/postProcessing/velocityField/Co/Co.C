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
    Co

Description
    Calculates and writes the Co number as a volScalarField obtained
    from field phi.

    The -noWrite option just outputs the max values without writing the
    field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject phiHeader
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (phiHeader.headerOk())
    {
        volScalarField Co
        (
            IOobject
            (
                "Co",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("0", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        Info<< "    Reading phi" << endl;
        surfaceScalarField phi(phiHeader, mesh);

        if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
        {
            Info<< "    Calculating compressible Co" << endl;

            Info<< "    Reading rho" << endl;
            volScalarField rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ
                ),
                mesh
            );

            Co.dimensionedInternalField() =
                (0.5*runTime.deltaT())
               *fvc::surfaceSum(mag(phi))().dimensionedInternalField()
               /(rho*mesh.V());
            Co.correctBoundaryConditions();
        }
        else if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
        {
            Info<< "    Calculating incompressible Co" << endl;

            Co.dimensionedInternalField() =
                (0.5*runTime.deltaT())
               *fvc::surfaceSum(mag(phi))().dimensionedInternalField()
               /mesh.V();
            Co.correctBoundaryConditions();
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Incorrect dimensions of phi: " << phi.dimensions()
                << abort(FatalError);
        }

        Info<< "Co max : " << max(Co).value() << endl;

        if (writeResults)
        {
            Co.write();
        }
    }
    else
    {
        Info<< "    No phi" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
