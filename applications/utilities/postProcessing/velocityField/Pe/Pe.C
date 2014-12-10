/*---------------------------------------------------------------------------* \
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
    Pe

Description
    Calculates and writes the Pe number as a surfaceScalarField obtained from
    field phi.

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"
#include "fluidThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "compressible/LES/LESModel/LESModel.H"

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
        autoPtr<surfaceScalarField> PePtr;

        Info<< "    Reading phi" << endl;
        surfaceScalarField phi(phiHeader, mesh);

        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        IOobject RASPropertiesHeader
        (
            "RASProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        );

        IOobject LESPropertiesHeader
        (
            "LESProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        );

        Info<< "    Calculating Pe" << endl;

        if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
        {
            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::RASModel> RASModel
                (
                    incompressible::RASModel::New
                    (
                        U,
                        phi,
                        laminarTransport
                    )
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(RASModel->nuEff())
                        )
                    )
                );
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                singlePhaseTransportModel laminarTransport(U, phi);

                autoPtr<incompressible::LESModel> sgsModel
                (
                    incompressible::LESModel::New(U, phi, laminarTransport)
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(sgsModel->nuEff())
                        )
                    )
                );
            }
            else
            {
                IOdictionary transportProperties
                (
                    IOobject
                    (
                        "transportProperties",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                dimensionedScalar nu(transportProperties.lookup("nu"));

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * nu
                        )
                    )
                );
            }
        }
        else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
        {
            if (RASPropertiesHeader.headerOk())
            {
                IOdictionary RASProperties(RASPropertiesHeader);

                autoPtr<fluidThermo> thermo(fluidThermo::New(mesh));

                volScalarField rho
                (
                    IOobject
                    (
                        "rho",
                        runTime.timeName(),
                        mesh
                    ),
                    thermo->rho()
                );

                autoPtr<compressible::RASModel> RASModel
                (
                    compressible::RASModel::New
                    (
                        rho,
                        U,
                        phi,
                        thermo()
                    )
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(RASModel->muEff())
                        )
                    )
                );
            }
            else if (LESPropertiesHeader.headerOk())
            {
                IOdictionary LESProperties(LESPropertiesHeader);

                autoPtr<fluidThermo> thermo(fluidThermo::New(mesh));

                volScalarField rho
                (
                    IOobject
                    (
                        "rho",
                        runTime.timeName(),
                        mesh
                    ),
                    thermo->rho()
                );

                autoPtr<compressible::LESModel> sgsModel
                (
                    compressible::LESModel::New(rho, U, phi, thermo())
                );

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * fvc::interpolate(sgsModel->muEff())
                        )
                    )
                );
            }
            else
            {
                IOdictionary transportProperties
                (
                    IOobject
                    (
                        "transportProperties",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                dimensionedScalar mu(transportProperties.lookup("mu"));

                PePtr.set
                (
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            "Pe",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        mag(phi)
                       /(
                            mesh.magSf()
                          * mesh.surfaceInterpolation::deltaCoeffs()
                          * mu
                        )
                    )
                );
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Incorrect dimensions of phi: " << phi.dimensions()
                    << abort(FatalError);
        }

        Info<< "Pe max : " << max(PePtr()).value() << endl;

        if (writeResults)
        {
            PePtr().write();
        }
    }
    else
    {
        Info<< "    No phi" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
