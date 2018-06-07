/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "drippingInjection.H"
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

defineTypeNameAndDebug(drippingInjection, 0);
addToRunTimeSelectionTable(injectionModel, drippingInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

drippingInjection::drippingInjection
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    deltaStable_(readScalar(coeffDict_.lookup("deltaStable"))),
    particlesPerParcel_(readScalar(coeffDict_.lookup("particlesPerParcel"))),
    rndGen_(label(0)),
    parcelDistribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("parcelDistribution"),
            rndGen_
        )
    ),
    diameter_(film.regionMesh().nCells(), -1.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

drippingInjection::~drippingInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void drippingInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const kinematicSingleLayer& film =
        refCast<const kinematicSingleLayer>(this->film());

    const scalar pi = constant::mathematical::pi;

    // calculate available dripping mass
    tmp<volScalarField> tgNorm(film.gNorm());
    const scalarField& gNorm = tgNorm();
    const scalarField& magSf = film.magSf();

    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();

    scalarField massDrip(film.regionMesh().nCells(), 0.0);

    forAll(gNorm, i)
    {
        if (gNorm[i] > small)
        {
            const scalar ddelta = max(0.0, delta[i] - deltaStable_);
            massDrip[i] +=
                min(availableMass[i], max(0.0, ddelta*rho[i]*magSf[i]));
        }
    }


    // Collect the data to be transferred
    forAll(massDrip, celli)
    {
        if (massDrip[celli] > 0)
        {
            // set new particle diameter if not already set
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
                massToInject[celli] = 0.0;
                diameterToInject[celli] = 0.0;
            }
        }
        else
        {
            massToInject[celli] = 0.0;
            diameterToInject[celli] = 0.0;
        }
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
