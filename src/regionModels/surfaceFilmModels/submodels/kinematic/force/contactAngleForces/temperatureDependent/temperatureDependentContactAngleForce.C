/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

#include "temperatureDependentContactAngleForce.H"
#include "thermoSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(temperatureDependentContactAngleForce, 0);
addToRunTimeSelectionTable
(
    force,
    temperatureDependentContactAngleForce,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temperatureDependentContactAngleForce::temperatureDependentContactAngleForce
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    contactAngleForce(typeName, film, dict),
    thetaPtr_(Function1<scalar>::New("theta", coeffDict_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temperatureDependentContactAngleForce::~temperatureDependentContactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> temperatureDependentContactAngleForce::theta() const
{
    tmp<volScalarField> ttheta
    (
        volScalarField::New
        (
            IOobject::modelName("theta", typeName),
            filmModel_.regionMesh(),
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField& theta = ttheta.ref();

    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    const volScalarField& T = film.T();

    theta.ref().field() = thetaPtr_->value(T());

    forAll(theta.boundaryField(), patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            theta.boundaryFieldRef()[patchi] =
                thetaPtr_->value(T.boundaryField()[patchi]);
        }
    }

    return ttheta;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
