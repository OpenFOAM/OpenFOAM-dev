/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Mach

Description
    Calculates and optionally writes the local Mach number from the velocity
    field U at each time.

    The -nowrite option just outputs the max value without writing the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fluidThermo.H"

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

    IOobject Theader
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check U and T exists
    if (Uheader.headerOk() && Theader.headerOk())
    {
        autoPtr<volScalarField> MachPtr;

        volVectorField U(Uheader, mesh);

        if
        (
            IOobject
            (
                "thermophysicalProperties",
                runTime.constant(),
                mesh
            ).headerOk()
        )
        {
            // thermophysical Mach
            autoPtr<fluidThermo> thermo
            (
                fluidThermo::New(mesh)
            );

            volScalarField Cp(thermo->Cp());
            volScalarField Cv(thermo->Cv());

            MachPtr.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "Ma",
                        runTime.timeName(),
                        mesh
                    ),
                    mag(U)/(sqrt((Cp/Cv)*(Cp - Cv)*thermo->T()))
                )
            );
        }
        else
        {
            // thermodynamic Mach
            IOdictionary thermoProps
            (
                IOobject
                (
                    "thermodynamicProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            );

            dimensionedScalar R(thermoProps.lookup("R"));
            dimensionedScalar Cv(thermoProps.lookup("Cv"));

            volScalarField T(Theader, mesh);

            MachPtr.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "Ma",
                        runTime.timeName(),
                        mesh
                    ),
                    mag(U)/(sqrt(((Cv + R)/Cv)*R*T))
                )
            );
        }

        Info<< "Mach max : " << max(MachPtr()).value() << endl;

        if (writeResults)
        {
            MachPtr().write();
        }
    }
    else
    {
        Info<< "    Missing U or T" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
