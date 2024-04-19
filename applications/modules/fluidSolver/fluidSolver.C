/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "fluidSolver.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(fluidSolver, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::fluidSolver::dependenciesModified() const
{
    return runTime.controlDict().modified() || mesh.solution().modified();
}


bool Foam::solvers::fluidSolver::read()
{
    solver::read();

    maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", vGreat);

    maxDeltaT_ =
        runTime.controlDict().found("maxDeltaT")
      ? runTime.controlDict().lookup<scalar>("maxDeltaT", runTime.userUnits())
      : vGreat;

    correctPhi = pimple.dict().lookupOrDefault
    (
        "correctPhi",
        mesh.dynamic()
    );

    checkMeshCourantNo = pimple.dict().lookupOrDefault
    (
        "checkMeshCourantNo",
        false
    );

    return true;
}

void Foam::solvers::fluidSolver::meshCourantNo() const
{
    if (checkMeshCourantNo)
    {
        const scalarField sumPhi
        (
            fvc::surfaceSum(mag(mesh.phi()))().primitiveField()
        );

        const scalar meshCoNum
        (
            0.5*gMax(sumPhi/mesh.V().primitiveField())*runTime.deltaTValue()
        );

        const scalar meanMeshCoNum
        (
            0.5
           *(gSum(sumPhi)/gSum(mesh.V().primitiveField()))
           *runTime.deltaTValue()
        );

        Info<< "Mesh Courant Number mean: " << meanMeshCoNum
            << " max: " << meshCoNum << endl;
    }
}


template<class RhoType>
void Foam::solvers::fluidSolver::correctCoNum
(
    const RhoType& rho,
    const surfaceScalarField& phi
)
{
    const scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()/rho.primitiveField()
    );

    CoNum_ = 0.5*gMax(sumPhi/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5
       *(gSum(sumPhi)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


void Foam::solvers::fluidSolver::correctCoNum(const surfaceScalarField& phi)
{
    correctCoNum(geometricOneField(), phi);
}


void Foam::solvers::fluidSolver::correctCoNum
(
    const volScalarField& rho,
    const surfaceScalarField& phi
)
{
    correctCoNum<volScalarField>(rho, phi);
}


void Foam::solvers::fluidSolver::continuityErrors
(
    const surfaceScalarField& phi
)
{
    const volScalarField contErr(fvc::div(phi));

    const scalar sumLocalContErr =
        runTime.deltaTValue()
       *mag(contErr)().weightedAverage(mesh.V()).value();

    const scalar globalContErr =
        runTime.deltaTValue()
       *contErr.weightedAverage(mesh.V()).value();

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr;

    if (pimple.finalPisoIter() && pimple.finalIter())
    {
        cumulativeContErr += globalContErr;

        Info<< ", cumulative = " << cumulativeContErr;
    }

    Info<< endl;
}


void Foam::solvers::fluidSolver::continuityErrors
(
    const volScalarField& rho,
    const volScalarField& thermoRho,
    const surfaceScalarField& phi
)
{
    if (mesh.schemes().steady())
    {
        continuityErrors(phi);
    }
    else
    {
        const dimensionedScalar totalMass = fvc::domainIntegrate(rho);

        const scalar sumLocalContErr =
            (fvc::domainIntegrate(mag(rho - thermoRho))/totalMass).value();

        const scalar globalContErr =
            (fvc::domainIntegrate(rho - thermoRho)/totalMass).value();

        cumulativeContErr += globalContErr;

        Info<< "time step continuity errors : sum local = " << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::fluidSolver::fluidSolver(fvMesh& mesh)
:
    solver(mesh),
    maxCo(0),
    maxDeltaT_(0),
    cumulativeContErr(0),
    CoNum_(0),
    CoNum(CoNum_)
{
    // Read the controls
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::fluidSolver::~fluidSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::fluidSolver::maxDeltaT() const
{
    scalar deltaT = min(fvModels().maxDeltaT(), maxDeltaT_);

    if (maxCo < vGreat && CoNum > small)
    {
        deltaT = min(deltaT, maxCo/CoNum*runTime.deltaTValue());
    }

    return deltaT;
}


// ************************************************************************* //
