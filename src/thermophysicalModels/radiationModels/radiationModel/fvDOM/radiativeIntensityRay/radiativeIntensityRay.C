/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "radiativeIntensityRay.H"
#include "fvm.H"
#include "fvDOM.H"
#include "constants.H"

using namespace Foam::constant;


const Foam::word
Foam::radiation::radiativeIntensityRay::intensityPrefix("ILambda");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::radiativeIntensityRay
(
    const fvDOM& dom,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const absorptionEmissionModel& absorptionEmission,
    const blackBodyEmission& blackBody,
    const label rayId
)
:
    dom_(dom),
    mesh_(mesh),
    absorptionEmission_(absorptionEmission),
    blackBody_(blackBody),
    I_
    (
        IOobject
        (
            "I" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("I", dimMass/pow3(dimTime), 0.0)
    ),
    Qr_
    (
        IOobject
        (
            "Qr" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    Qin_
    (
        IOobject
        (
            "Qin" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimMass/pow3(dimTime), 0.0)
    ),
    Qem_
    (
        IOobject
        (
            "Qem" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qem", dimMass/pow3(dimTime), 0.0)
    ),
    d_(vector::zero),
    dAve_(vector::zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    nLambda_(nLambda),
    ILambda_(nLambda),
    myRayId_(rayId)
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    d_ = vector(sinTheta*sinPhi, sinTheta*cosPhi, cosTheta);
    dAve_ = vector
    (
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        cosPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        0.5*deltaPhi*Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
    );

    autoPtr<volScalarField> IDefaultPtr;

    forAll(ILambda_, lambdaI)
    {
        IOobject IHeader
        (
            intensityPrefix + "_" + name(rayId) + "_" + name(lambdaI),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // check if field exists and can be read
        if (IHeader.headerOk())
        {
            ILambda_.set
            (
                lambdaI,
                new volScalarField(IHeader, mesh_)
            );
        }
        else
        {
            // Demand driven load the IDefault field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "IDefault",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }

            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

            ILambda_.set
            (
                lambdaI,
                new volScalarField(noReadHeader, IDefaultPtr())
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::radiativeIntensityRay::~radiativeIntensityRay()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::radiativeIntensityRay::correct()
{
    // reset boundary heat flux to zero
    Qr_.boundaryField() = 0.0;

    scalar maxResidual = -GREAT;

    forAll(ILambda_, lambdaI)
    {
        const volScalarField& k = dom_.aLambda(lambdaI);

        tmp<fvScalarMatrix> IiEq;

        if (!dom_.cacheDiv())
        {
            const surfaceScalarField Ji(dAve_ & mesh_.Sf());

            IiEq =
            (
                fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
              + fvm::Sp(k*omega_, ILambda_[lambdaI])
            ==
                1.0/constant::mathematical::pi*omega_
              * (
                    k*blackBody_.bLambda(lambdaI)
                  + absorptionEmission_.ECont(lambdaI)/4
                )
            );
        }
        else
        {
            IiEq =
            (
               dom_.fvRayDiv(myRayId_, lambdaI)
             + fvm::Sp(k*omega_, ILambda_[lambdaI])
           ==
               1.0/constant::mathematical::pi*omega_
             * (
                   k*blackBody_.bLambda(lambdaI)
                 + absorptionEmission_.ECont(lambdaI)/4
               )
            );
        }

        IiEq().relax();

        const solverPerformance ILambdaSol = solve
        (
            IiEq(),
            mesh_.solver("Ii")
        );

        const scalar initialRes =
            ILambdaSol.initialResidual()*omega_/dom_.omegaMax();

        maxResidual = max(initialRes, maxResidual);
    }

    return maxResidual;
}


void Foam::radiation::radiativeIntensityRay::addIntensity()
{
    I_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);

    forAll(ILambda_, lambdaI)
    {
        I_ += ILambda_[lambdaI];
    }
}


// ************************************************************************* //
