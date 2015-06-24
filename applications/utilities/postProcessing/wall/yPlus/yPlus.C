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
    yPlus

Description
    Calculates and reports yPlus for the near-wall cells of all wall patches,
    for the specified times for laminar, LES and RAS.

    For walls at which wall-functions are applied the wall-function provides
    the y+ values otherwise they are obtained directly from the near-wall
    velocity gradient and effective and laminar viscosities.

    Compressible modes is automatically selected based on the existence of the
    "thermophysicalProperties" dictionary required to construct the
    thermodynamics package.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TurbulenceModel>
void calcYPlus
(
    const TurbulenceModel& turbulenceModel,
    const fvMesh& mesh,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();

    const volScalarField::GeometricBoundaryField nutBf =
        turbulenceModel->nut()().boundaryField();

    const volScalarField::GeometricBoundaryField nuEffBf =
        turbulenceModel->nuEff()().boundaryField();

    const volScalarField::GeometricBoundaryField nuBf =
        turbulenceModel->nu()().boundaryField();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<nutWallFunctionFvPatchScalarField>(nutBf[patchi]))
        {
            const nutWallFunctionFvPatchScalarField& nutPf =
                dynamic_cast<const nutWallFunctionFvPatchScalarField&>
                (
                    nutBf[patchi]
                );

            yPlus.boundaryField()[patchi] = nutPf.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << nutPf.patch().name()
                << ", wall-function " << nutPf.type()
                << ", y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
        else if (isA<wallFvPatch>(patch))
        {
            yPlus.boundaryField()[patchi] =
                d[patchi]
               *sqrt
                (
                    nuEffBf[patchi]
                   *mag(U.boundaryField()[patchi].snGrad())
                )/nuBf[patchi];
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << patch.name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }
}


void calcIncompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulenceModel
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    calcYPlus(turbulenceModel, mesh, U, yPlus);
}


void calcCompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
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

    autoPtr<compressible::turbulenceModel> turbulenceModel
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    calcYPlus(turbulenceModel, mesh, U, yPlus);
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select(runTime, args, "yPlus");
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
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
                calcCompressibleYPlus(mesh, runTime, U, yPlus);
            }
            else
            {
                calcIncompressibleYPlus(mesh, runTime, U, yPlus);
            }
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "Writing yPlus to field " << yPlus.name() << nl << endl;

        yPlus.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
