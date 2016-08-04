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
    dCoeff_(coeffDict_.lookupOrDefault("dCoeff", 3.3)),
    deltaStable_(coeffDict_.lookupOrDefault("deltaStable", 0)),
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

    // Calculate available dripping mass
    tmp<volScalarField> tsinAlpha(film.gNorm()/mag(film.g()));
    const scalarField& sinAlpha = tsinAlpha();
    const scalarField& magSf = film.magSf();

    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();
    const scalarField& sigma = film.sigma();
    const scalar magg = mag(film.g().value());

    scalarField massDrip(film.regionMesh().nCells(), scalar(0));

    forAll(delta, celli)
    {
        if (sinAlpha[celli] > SMALL && delta[celli] > deltaStable_)
        {
            const scalar lc = sqrt(sigma[celli]/(rho[celli]*magg));
            const scalar deltaStable = max
            (
                3*lc*sqrt(1 - sqr(sinAlpha[celli]))
               /(ubarStar_*sqrt(sinAlpha[celli])*sinAlpha[celli]),
                deltaStable_
            );

            if (delta[celli] > deltaStable)
            {
                const scalar ddelta = max(delta[celli] - deltaStable, 0);
                massDrip[celli] +=
                    min
                    (
                        availableMass[celli],
                        max(ddelta*rho[celli]*magSf[celli], 0)
                    );
            }
        }
    }

    // Collect the data to be transferred
    forAll(massDrip, celli)
    {
        if (massDrip[celli] > 0)
        {
            const scalar rhoc = rho[celli];
            const scalar diam = dCoeff_*sqrt(sigma[celli]/(rhoc*magg));
            diameter_[celli] = diam;

            massToInject[celli] += massDrip[celli];
            availableMass[celli] -= massDrip[celli];

            diameterToInject[celli] = diam;
            addToInjectedMass(massDrip[celli]);
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
