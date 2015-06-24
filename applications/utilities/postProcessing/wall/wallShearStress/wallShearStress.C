/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    wallShearStress

Description
    Calculates and reports wall shear stress for all patches, for the
    specified times when using RAS turbulence models.

    Compressible modes is automatically selected based on the existence of the
    "thermophysicalProperties" dictionary required to construct the
    thermodynamics package.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& wallShearStress
)
{
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> model
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    const volSymmTensorField Reff(model->devReff());

    forAll(wallShearStress.boundaryField(), patchI)
    {
        wallShearStress.boundaryField()[patchI] =
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI];
    }
}


void calcCompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volVectorField& wallShearStress
)
{
    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!rhoHeader.headerOk())
    {
        Info<< "    no rho field" << endl;
        return;
    }

    Info<< "Reading field rho\n" << endl;
    volScalarField rho(rhoHeader, mesh);

    #include "compressibleCreatePhi.H"

    autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
    fluidThermo& thermo = pThermo();

    autoPtr<compressible::RASModel> model
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    const volSymmTensorField Reff(model->devRhoReff());

    forAll(wallShearStress.boundaryField(), patchI)
    {
        wallShearStress.boundaryField()[patchI] =
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI];
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volVectorField wallShearStress
        (
            IOobject
            (
                "wallShearStress",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "wallShearStress",
                sqr(dimLength)/sqr(dimTime),
                vector::zero
            )
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);

            if
            (
                IOobject
                (
                    basicThermo::dictName,
                    runTime.constant(),
                    mesh
                ).headerOk()
            )
            {
                calcCompressible(mesh, runTime, U, wallShearStress);
            }
            else
            {
                calcIncompressible(mesh, runTime, U, wallShearStress);
            }
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "Writing wall shear stress to field " << wallShearStress.name()
            << nl << endl;

        wallShearStress.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
