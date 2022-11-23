/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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

#include "perturbedTemperatureDependentContactAngleForce.H"
#include "thermoSurfaceFilm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(perturbedTemperatureDependentContactAngleForce, 0);
addToRunTimeSelectionTable
(
    force,
    perturbedTemperatureDependentContactAngleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

perturbedTemperatureDependentContactAngleForce::
perturbedTemperatureDependentContactAngleForce
(
    surfaceFilm& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    thetaPtr_(Function1<scalar>::New("theta", coeffDict_)),
    rndGen_(label(0)),
    distribution_
    (
        distributionModel::New
        (
            coeffDict_.subDict("distribution"),
            rndGen_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

perturbedTemperatureDependentContactAngleForce::
~perturbedTemperatureDependentContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField>
perturbedTemperatureDependentContactAngleForce::theta() const
{
    tmp<volScalarField> ttheta
    (
        volScalarField::New
        (
            typedName("theta"),
            filmModel_.mesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField& theta = ttheta.ref();
    volScalarField::Internal& thetai = theta.ref();

    const thermoSurfaceFilm& film = filmType<thermoSurfaceFilm>();

    const volScalarField& T = film.thermo().T();

    // Initialise with the function of temperature
    thetai.field() = thetaPtr_->value(T());

    // Add the stochastic perturbation
    forAll(thetai, celli)
    {
        thetai[celli] += distribution_->sample();
    }

    forAll(theta.boundaryField(), patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            fvPatchField<scalar>& thetaf = theta.boundaryFieldRef()[patchi];

            // Initialise with the function of temperature
            thetaf = thetaPtr_->value(T.boundaryField()[patchi]);

            // Add the stochastic perturbation
            forAll(thetaf, facei)
            {
                thetaf[facei] += distribution_->sample();
            }
        }
    }

    return ttheta;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// ************************************************************************* //
