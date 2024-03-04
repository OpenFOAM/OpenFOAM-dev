/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "BrunDripping.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmEjectionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(BrunDripping, 0);
addToRunTimeSelectionTable(ejectionModel, BrunDripping, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BrunDripping::BrunDripping
(
    const dictionary& dict,
    const solvers::isothermalFilm& film
)
:
    ejectionModel(dict, film),
    ubarStar_
    (
        dict.optionalSubDict(typeName + "Coeffs")
       .lookupOrDefault("ubarStar", 1.62208)
    ),
    dCoeff_
    (
        dict.optionalSubDict(typeName + "Coeffs")
       .lookupOrDefault("dCoeff", 3.3)
    ),
    deltaStable_
    (
        dict.optionalSubDict(typeName + "Coeffs")
       .lookupOrDefault("deltaStable", scalar(0))
    ),
    minParticlesPerParcel_
    (
        dict.optionalSubDict(typeName + "Coeffs")
       .lookupOrDefault("minParticlesPerParcel", 1)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

BrunDripping::~BrunDripping()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void BrunDripping::correct()
{
    const scalar piBy6 = constant::mathematical::pi/6;

    const scalarField& magSf = film_.magSf;
    const vectorField& nHat = film_.nHat;
    const scalarField& delta = film_.delta;
    const scalarField& rho = film_.rho;

    const tmp<volScalarField> tsigma = film_.sigma();
    const volScalarField::Internal& sigma = tsigma();

    const scalar magg = mag(film_.g.value());
    const vector gHat = -film_.g.value()/magg;

    const scalar deltaT = film_.mesh.time().deltaTValue();

    forAll(delta, celli)
    {
        rate_[celli] = 0;
        diameter_[celli] = 0;

        const scalar sinAlpha = nHat[celli] & gHat;

        if (sinAlpha > small && delta[celli] > deltaStable_)
        {
            const scalar lc = sqrt(sigma[celli]/(rho[celli]*magg));

            const scalar deltaStable = max
            (
                3*lc*sqrt(1 - sqr(sinAlpha))
               /(ubarStar_*sqrt(sinAlpha)*sinAlpha),
                deltaStable_
            );

            if (delta[celli] > deltaStable)
            {
                const scalar ddelta = delta[celli] - deltaStable;
                const scalar massDrip = ddelta*rho[celli]*magSf[celli];

                // Calculate dripped droplet diameter
                diameter_[celli] = dCoeff_*lc;

                // Calculate the minimum mass of a droplet parcel
                const scalar minMass =
                    minParticlesPerParcel_
                   *rho[celli]*piBy6*pow3(diameter_[celli]);

                if (massDrip > minMass)
                {
                    // Calculate rate of dripping
                    rate_[celli] = ddelta/(deltaT*delta[celli]);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmEjectionModels
} // End namespace Foam

// ************************************************************************* //
