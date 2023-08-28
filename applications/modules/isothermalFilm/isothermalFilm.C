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
#include "mappedPatchBase.H"
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
            << "The number of film wall faces in the mesh "
            << nWallFaces
            << " is not equal to the number of cells "
            << mesh.nCells()
            << exit(FatalError);
    }

    if (returnReduce(nWallFaces, sumOp<label>()) == 0)
    {
        FatalErrorInFunction
            << "There are no filmWall faces in the mesh"
            << exit(FatalError);
    }

    wallPatchIDs.transfer(wallPatchIDs_);


    // Search for film surface patch
    surfacePatchID = -1;

    forAll(bm, patchi)
    {
        const polyPatch& p = bm[patchi];

        if (isA<filmSurfacePolyPatch>(p))
        {
            if (surfacePatchID == -1)
            {
                surfacePatchID = patchi;
            }
            else
            {
                FatalErrorInFunction
                    << "More than one filmSurface patch defined: "
                    << surfacePatchID << " and " << patchi
                    << exit(FatalError);
            }
        }
    }

    if (surfacePatchID == -1)
    {
        Info<< "The filmSurface patch is not defined"
            << endl;
    }


    // Calculate film specific mesh geometry from the film wall patches

    forAll(wallPatchIDs, i)
    {
        const label patchi = wallPatchIDs[i];
        const polyPatch& wallp = bm[patchi];
        const labelList& fCells = wallp.faceCells();

        UIndirectList<vector>(nHat_, fCells) = wallp.faceNormals();
        UIndirectList<scalar>(magSf_, fCells) = wallp.magFaceAreas();
    }

    nHat_.correctBoundaryConditions();

    VbyA_.primitiveFieldRef() = mesh.V()/magSf_;
    VbyA_.correctBoundaryConditions();

    return true;
}


Foam::wordList Foam::solvers::isothermalFilm::alphaTypes() const
{
    wordList alphaTypes(delta_.boundaryField().types());

    forAll(delta_.boundaryField(), patchi)
    {
        if (!delta_.boundaryField()[patchi].assignable())
        {
            alphaTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    forAll(wallPatchIDs, i)
    {
        alphaTypes[wallPatchIDs[i]] = alphaOneFvPatchScalarField::typeName;
    }

    if (surfacePatchID != -1)
    {
        alphaTypes[surfacePatchID] = alphaOneFvPatchScalarField::typeName;
    }

    return alphaTypes;
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

    correctContinuityError();

    if (mass.value() > small)
    {
        const volScalarField::Internal massContErr
        (
            runTime.deltaT()*magSf*contErr()
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

bool Foam::solvers::isothermalFilm::dependenciesModified() const
{
    return runTime.controlDict().modified();
}


bool Foam::solvers::isothermalFilm::read()
{
    solver::read();

    maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", vGreat);

    maxDeltaT_ =
        runTime.controlDict().found("maxDeltaT")
      ? runTime.userTimeToTime
        (
            runTime.controlDict().lookup<scalar>("maxDeltaT")
        )
      : vGreat;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::isothermalFilm::isothermalFilm
(
    fvMesh& mesh,
    autoPtr<rhoFluidThermo> thermoPtr
)
:
    solver(mesh),

    CoNum(0),
    cumulativeContErr(0),

    thermoPtr_(thermoPtr),
    thermo_(thermoPtr_()),

    p(thermo_.p()),

    nHat_
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

    magSf_
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

    VbyA_
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

    delta_
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

    alpha_
    (
        IOobject
        (
            "alpha",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        delta_/VbyA_,
        alphaTypes()
    ),

    deltaWet("deltaWet", dimLength, thermo_.properties()),

    U_
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

    alphaRhoPhi_
    (
        IOobject
        (
            "alphaRhoPhi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(alpha_*thermo_.rho()*U_)
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U_)
    ),

    surfaceTension(surfaceTensionModel::New(thermo_.properties(), mesh)),
    thermocapillary(!isType<surfaceTensionModels::constant>(surfaceTension())),

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

    nHat(nHat_),
    magSf(magSf_),
    VbyA(VbyA_),
    delta(delta_),
    alpha(alpha_),
    thermo(thermo_),
    rho(thermo_.rho()),
    U(U_),
    alphaRhoPhi(alphaRhoPhi_),
    phi(phi_),

    momentumTransport
    (
        filmCompressible::momentumTransportModel::New
        (
            alpha,
            thermo.rho(),
            U,
            alphaRhoPhi,
            phi,
            thermo
        )
    )
{
    // Read the controls
    read();

    mesh.schemes().setFluxRequired(alpha.name());
    momentumTransport->validate();
}


Foam::solvers::isothermalFilm::isothermalFilm(fvMesh& mesh)
:
    isothermalFilm(mesh, rhoFluidThermo::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::isothermalFilm::~isothermalFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::fvPatch& Foam::solvers::isothermalFilm::surfacePatch() const
{
    return mesh.boundary()[surfacePatchID];
}


const Foam::mappedPatchBase&
Foam::solvers::isothermalFilm::surfacePatchMap() const
{
    return refCast<const mappedPatchBase>(surfacePatch().patch());
}


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
    correctCoNum();
}


void Foam::solvers::isothermalFilm::moveMesh()
{}


void Foam::solvers::isothermalFilm::thermophysicalPredictor()
{
    thermo_.correct();
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
