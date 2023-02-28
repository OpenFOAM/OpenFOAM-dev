/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "isothermalFilm.H"
#include "filmWallPolyPatch.H"
#include "filmSurfacePolyPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "alphaOneFvPatchScalarField.H"
#include "constantSurfaceTension.H"
#include "fvcVolumeIntegrate.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcSurfaceIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(isothermalFilm, 0);
    addToRunTimeSelectionTable(solver, isothermalFilm, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::isothermalFilm::readControls()
{
    maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxDeltaT_ =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", vGreat);
}


void Foam::solvers::isothermalFilm::correctCoNum()
{
    const scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


void Foam::solvers::isothermalFilm::continuityErrors()
{
    const dimensionedScalar mass = fvc::domainIntegrate(rho()*delta()*magSf);

    if (mass.value() > small)
    {
        const volScalarField::Internal contErr
        (
            fvc::ddt(alpha, rho)()() + fvc::div(alphaRhoPhi)()()
          - (fvModels().source(rho, alpha) & alpha)()()
        );

        const volScalarField::Internal massContErr
        (
            runTime.deltaT()*magSf*contErr
        );

        const scalar sumLocalContErr =
            (fvc::domainIntegrate(mag(massContErr))/mass).value();

        const scalar globalContErr =
            (fvc::domainIntegrate(massContErr)/mass).value();

        Info<< "time step continuity errors : sum local = " << sumLocalContErr
            << ", global = " << globalContErr;

        if (pimple.finalPisoIter() && pimple.finalIter())
        {
            cumulativeContErr += globalContErr;

            Info<< ", cumulative = " << cumulativeContErr;
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::isothermalFilm::initFilmMesh()
{
    // Search for film wall patches

    label nWallFaces = 0;
    DynamicList<label> wallPatchIDs_;

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    forAll(bm, patchi)
    {
        const polyPatch& p = bm[patchi];

        if (isA<filmWallPolyPatch>(p))
        {
            wallPatchIDs_.append(patchi);
            nWallFaces += p.faceCells().size();
        }
    }

    if (nWallFaces != mesh.nCells())
    {
        FatalErrorInFunction
            << "The number of filmWall faces in the mesh "
               "is not equal to the number of cells"
            << exit(FatalError);
    }

    if (returnReduce(nWallFaces, sumOp<label>()) == 0)
    {
        FatalErrorInFunction
            << "There are no filmWall faces in the mesh"
            << exit(FatalError);
    }

    wallPatchIDs.transfer(wallPatchIDs_);


    // Search for film surface patches

    label nSurfaceFaces = 0;
    DynamicList<label> surfacePatchIDs_;

    forAll(bm, patchi)
    {
        const polyPatch& p = bm[patchi];

        if (isA<filmSurfacePolyPatch>(p))
        {
            surfacePatchIDs_.append(patchi);
            nSurfaceFaces += p.faceCells().size();
        }
    }

    surfacePatchIDs.transfer(surfacePatchIDs_);

    if (returnReduce(nSurfaceFaces, sumOp<label>()) == 0)
    {
        Info<< "There are no filmSurface faces in the mesh"
            << endl;
    }


    // Calculate film specific mesh geometry from the film wall patches

    forAll(wallPatchIDs, i)
    {
        const label patchi = wallPatchIDs[i];
        const polyPatch& wallp = bm[patchi];
        const labelList& fCells = wallp.faceCells();

        UIndirectList<vector>(nHat, fCells) = wallp.faceNormals();
        UIndirectList<scalar>(magSf, fCells) = wallp.magFaceAreas();
    }

    nHat.correctBoundaryConditions();

    VbyA.primitiveFieldRef() = mesh.V()/magSf;
    VbyA.correctBoundaryConditions();

    return true;
}


Foam::wordList Foam::solvers::isothermalFilm::alphaTypes() const
{
    wordList alphaTypes(delta.boundaryField().types());

    forAll(delta.boundaryField(), patchi)
    {
        if (!delta.boundaryField()[patchi].assignable())
        {
            alphaTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    forAll(wallPatchIDs, i)
    {
        alphaTypes[wallPatchIDs[i]] = alphaOneFvPatchScalarField::typeName;
    }

    forAll(surfacePatchIDs, i)
    {
        alphaTypes[surfacePatchIDs[i]] = alphaOneFvPatchScalarField::typeName;
    }

    return alphaTypes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::isothermalFilm::isothermalFilm
(
    fvMesh& mesh,
    autoPtr<fluidThermo> thermoPtr
)
:
    solver(mesh),

    CoNum(0),
    cumulativeContErr(0),

    thermo_(thermoPtr),
    thermo(thermo_()),

    rho(thermo.rho()),

    nHat
    (
        IOobject
        (
            "nHat",
            runTime.name(),
            mesh
        ),
        mesh,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchField<vector>::typeName
    ),

    magSf
    (
        IOobject
        (
            "magSf",
            runTime.name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimArea, 0)
    ),

    VbyA
    (
        IOobject
        (
            "VbyA",
            runTime.name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimLength, 0),
        zeroGradientFvPatchField<vector>::typeName
    ),

    initialised_(initFilmMesh()),

    delta
    (
        IOobject
        (
            "delta",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alpha
    (
        IOobject
        (
            "alpha",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        delta/VbyA,
        alphaTypes()
    ),

    deltaWet("deltaWet", dimLength, thermo.properties()),

    g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alphaRhoPhi
    (
        IOobject
        (
            "alphaRhoPhi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(alpha*rho*U)
    ),

    phi
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    ),

    surfaceTension(surfaceTensionModel::New(thermo.properties(), mesh)),
    thermocapillary(!isType<surfaceTensionModels::constant>(surfaceTension())),

    momentumTransport
    (
        filmCompressible::momentumTransportModel::New
        (
            alpha,
            rho,
            U,
            alphaRhoPhi,
            phi,
            thermo
        )
    )
{
    // Read the controls
    readControls();

    mesh.schemes().setFluxRequired(alpha.name());
    momentumTransport->validate();
}


Foam::solvers::isothermalFilm::isothermalFilm(fvMesh& mesh)
:
    isothermalFilm(mesh, fluidThermo::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::isothermalFilm::~isothermalFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::isothermalFilm::maxDeltaT() const
{
    scalar deltaT = min(fvModels().maxDeltaT(), maxDeltaT_);

    if (CoNum > small)
    {
        deltaT = min(deltaT, maxCo/CoNum*runTime.deltaTValue());
    }

    return deltaT;
}


void Foam::solvers::isothermalFilm::preSolve()
{
    readControls();
    correctCoNum();
}


void Foam::solvers::isothermalFilm::moveMesh()
{}


void Foam::solvers::isothermalFilm::thermophysicalPredictor()
{
    thermo.correct();
}


void Foam::solvers::isothermalFilm::pressureCorrector()
{
    correctAlpha();
}


void Foam::solvers::isothermalFilm::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport->correct();
    }
}


void Foam::solvers::isothermalFilm::postSolve()
{}


// ************************************************************************* //
