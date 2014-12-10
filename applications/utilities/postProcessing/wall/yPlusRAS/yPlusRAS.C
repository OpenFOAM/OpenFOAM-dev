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
    yPlusRAS

Description
    Calculates and reports yPlus for all wall patches, for the specified times
    when using RAS turbulence models.

    Default behaviour assumes operating in incompressible mode.
    Use the -compressible option for compressible RAS cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "nutWallFunction/nutWallFunctionFvPatchScalarField.H"

#include "fluidThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "mutWallFunction/mutWallFunctionFvPatchScalarField.H"

#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    typedef incompressible::nutWallFunctionFvPatchScalarField
        wallFunctionPatchField;

    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> RASModel
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    const volScalarField::GeometricBoundaryField nutPatches =
        RASModel->nut()().boundaryField();

    bool foundNutPatch = false;
    forAll(nutPatches, patchi)
    {
        if (isA<wallFunctionPatchField>(nutPatches[patchi]))
        {
            foundNutPatch = true;

            const wallFunctionPatchField& nutPw =
                dynamic_cast<const wallFunctionPatchField&>
                    (nutPatches[patchi]);

            yPlus.boundaryField()[patchi] = nutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << nutPw.patch().name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }

    if (!foundNutPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


void calcCompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    typedef compressible::mutWallFunctionFvPatchScalarField
        wallFunctionPatchField;

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

    autoPtr<fluidThermo> pThermo
    (
        fluidThermo::New(mesh)
    );
    fluidThermo& thermo = pThermo();

    autoPtr<compressible::RASModel> RASModel
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    const volScalarField::GeometricBoundaryField mutPatches =
        RASModel->mut()().boundaryField();

    bool foundMutPatch = false;
    forAll(mutPatches, patchi)
    {
        if (isA<wallFunctionPatchField>(mutPatches[patchi]))
        {
            foundMutPatch = true;

            const wallFunctionPatchField& mutPw =
                dynamic_cast<const wallFunctionPatchField&>
                    (mutPatches[patchi]);

            yPlus.boundaryField()[patchi] = mutPw.yPlus();
            const scalarField& Yp = yPlus.boundaryField()[patchi];

            Info<< "Patch " << patchi
                << " named " << mutPw.patch().name()
                << " y+ : min: " << gMin(Yp) << " max: " << gMax(Yp)
                << " average: " << gAverage(Yp) << nl << endl;
        }
    }

    if (!foundMutPatch)
    {
        Info<< "    no " << wallFunctionPatchField::typeName << " patches"
            << endl;
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "compressible",
        "calculate compressible y+"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const bool compressible = args.optionFound("compressible");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        // Wall distance
        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" << endl;
            wallDist y(mesh, true);
            Info<< "Writing wall distance to field " << y.name() << nl << endl;
            y.write();
        }

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

            if (compressible)
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
