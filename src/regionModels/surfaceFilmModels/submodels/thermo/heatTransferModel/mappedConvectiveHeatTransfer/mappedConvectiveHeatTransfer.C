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

#include "mappedConvectiveHeatTransfer.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumSurfaceFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmSubModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mappedConvectiveHeatTransfer, 0);

addToRunTimeSelectionTable
(
    heatTransferModel,
    mappedConvectiveHeatTransfer,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mappedConvectiveHeatTransfer::mappedConvectiveHeatTransfer
(
    surfaceFilm& film,
    const dictionary& dict
)
:
    heatTransferModel(film),
    htcConvPrimary_
    (
        IOobject
        (
            "htcConv",
            film.time().timeName(),
            film.primaryMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.primaryMesh()
    ),
    htcConvFilm_
    (
        IOobject
        (
            htcConvPrimary_.name(),
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime)/dimTemperature, 0)
    )
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mappedConvectiveHeatTransfer::~mappedConvectiveHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void mappedConvectiveHeatTransfer::correct()
{
    // Update the primary-side convective heat transfer coefficient
    htcConvPrimary_.correctBoundaryConditions();

    // Map the primary-side convective heat transfer coefficient
    // to the region internal field
    film().toRegion(htcConvFilm_, htcConvPrimary_.boundaryField());
}


tmp<volScalarField::Internal> mappedConvectiveHeatTransfer::h() const
{
    return htcConvFilm_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
