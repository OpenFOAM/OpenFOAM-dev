/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "dripping.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmEjectionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dripping, 0);
addToRunTimeSelectionTable(ejectionModel, dripping, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dripping::dripping
(
    const dictionary& dict,
    const solvers::isothermalFilm& film
)
:
    ejectionModel(dict, film),
    deltaStable_
    (
        dict.optionalSubDict(typeName + "Coeffs")
       .lookup<scalar>("deltaStable")
    ),
    minParticlesPerParcel_
    (
        dict.optionalSubDict(typeName + "Coeffs")
       .lookupOrDefault("minParticlesPerParcel", 1)
    ),
    rndGen_(label(0)),
    parcelDistribution_
    (
        distribution::New
        (
            dict.optionalSubDict(typeName + "Coeffs")
           .subDict("parcelDistribution"),
            rndGen_,
            0
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dripping::~dripping()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void dripping::correct()
{
    const scalar piBy6 = constant::mathematical::pi/6;

    const scalarField gnHat(film_.nHat() & (-film_.g));
    const scalarField& magSf = film_.magSf;

    const scalarField& delta = film_.delta;
    const scalarField& rho = film_.rho;

    const scalar deltaT = film_.mesh.time().deltaTValue();

    forAll(delta, celli)
    {
        rate_[celli] = 0;
        diameter_[celli] = 0;

        // Calculate available dripping mass
        if (gnHat[celli] > small && delta[celli] > deltaStable_)
        {
            const scalar ddelta = delta[celli] - deltaStable_;
            const scalar massDrip = ddelta*rho[celli]*magSf[celli];

            // Sample dripped droplet diameter
            diameter_[celli] = parcelDistribution_->sample();

            // Calculate the minimum mass of a parcel
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmEjectionModels
} // End namespace Foam

// ************************************************************************* //
