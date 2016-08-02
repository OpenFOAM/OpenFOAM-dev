/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "BrunDrippingInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "volFields.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BrunDrippingInjection, 0);
addToRunTimeSelectionTable(injectionModel, BrunDrippingInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BrunDrippingInjection::BrunDrippingInjection
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    injectionModel(type(), owner, dict),
    ubarStar_(coeffDict_.lookupOrDefault("ubarStar", 1.62208)),
    deltaStable_(readScalar(coeffDict_.lookup("deltaStable"))),
    particlesPerParcel_(readScalar(coeffDict_.lookup("particlesPerParcel"))),
    rndGen_(label(0), -1),
    parcelDistribution_
    (
        distributionModels::distributionModel::New
        (
            coeffDict_.subDict("parcelDistribution"),
            rndGen_
        )
    ),
    diameter_(owner.regionMesh().nCells(), -1.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

BrunDrippingInjection::~BrunDrippingInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void BrunDrippingInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const kinematicSingleLayer& film =
        refCast<const kinematicSingleLayer>(this->owner());

    const scalar pi = constant::mathematical::pi;

    // Calculate available dripping mass
    tmp<volScalarField> tsinAlpha(film.gNorm()/mag(film.g()));
    const scalarField& sinAlpha = tsinAlpha();
    const scalarField& magSf = film.magSf();

    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();
    const scalarField& sigma = film.sigma();
    const scalar magg = mag(film.g().value());

    scalarField massDrip(film.regionMesh().nCells(), 0.0);

    forAll(delta, i)
    {
        if (sinAlpha[i] > SMALL && delta[i] > deltaStable_)
        {
            scalar lc = sqrt(sigma[i]/(rho[i]*magg));
            scalar deltaStable = max
            (
                3*lc*sqrt(1 - sqr(sinAlpha[i]))
               /(ubarStar_*sqrt(sinAlpha[i])*sinAlpha[i]),
                deltaStable_
            );

            Info<< delta[i] << " " << deltaStable << endl;
            if (delta[i] > deltaStable)
            {
                const scalar ddelta = max(0.0, delta[i] - deltaStable);
                massDrip[i] +=
                    min(availableMass[i], max(0.0, ddelta*rho[i]*magSf[i]));
            }
        }
    }


    // Collect the data to be transferred
    forAll(massDrip, celli)
    {
        if (massDrip[celli] > 0)
        {
            // Set new particle diameter if not already set
            if (diameter_[celli] < 0)
            {
                diameter_[celli] = parcelDistribution_->sample();
            }

            scalar& diam = diameter_[celli];
            scalar rhoc = rho[celli];
            scalar minMass = particlesPerParcel_*rhoc*pi/6*pow3(diam);

            if (massDrip[celli] > minMass)
            {
                // All drip mass can be injected
                massToInject[celli] += massDrip[celli];
                availableMass[celli] -= massDrip[celli];

                // Set particle diameter
                diameterToInject[celli] = diam;

                // Retrieve new particle diameter sample
                diam = parcelDistribution_->sample();

                addToInjectedMass(massDrip[celli]);
            }
            else
            {
                // Particle mass below minimum threshold - cannot be injected
                massToInject[celli] = 0;
                diameterToInject[celli] = 0;
            }
        }
        else
        {
            massToInject[celli] = 0;
            diameterToInject[celli] = 0;
        }
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
