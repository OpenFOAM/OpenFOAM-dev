/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "standardRadiation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRadiation, 0);

addToRunTimeSelectionTable
(
    radiationModel,
    standardRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardRadiation::standardRadiation
(
     surfaceFilm& film,
    const dictionary& dict
)
:
    radiationModel(typeName, film, dict),
    qinFilm_
    (
        IOobject
        (
            "qin",
            film.time().name(),
            film.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.mesh(),
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    beta_(coeffDict_.lookup<scalar>("beta")),
    kappaBar_(coeffDict_.lookup<scalar>("kappaBar"))
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRadiation::~standardRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardRadiation::correct()
{
    const volScalarField& qinPrimary
    (
        film().primaryMesh().lookupObject<volScalarField>("qin")
    );

    // Map the primary-side radiative flux to the region internal field
    film().toRegion(qinFilm_, qinPrimary.boundaryField());
}


tmp<volScalarField::Internal> standardRadiation::Shs()
{
    return volScalarField::Internal::New
    (
        typedName("Shs"),
        beta_*qinFilm_*filmModel_.coverage()
       *(1 - exp(-kappaBar_*filmModel_.delta()()))
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// ************************************************************************* //
