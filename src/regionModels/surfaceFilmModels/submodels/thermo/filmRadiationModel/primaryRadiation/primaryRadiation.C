/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "primaryRadiation.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(primaryRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    primaryRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

primaryRadiation::primaryRadiation
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    filmRadiationModel(typeName, owner, dict),
    QinPrimary_
    (
        IOobject
        (
            "Qin", // same name as Qin on primary region to enable mapping
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0),
        owner.mappedPushedFieldPatchTypes<scalar>()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

primaryRadiation::~primaryRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void primaryRadiation::correct()
{
    // Transfer Qin from primary region
    QinPrimary_.correctBoundaryConditions();
}


tmp<volScalarField> primaryRadiation::Shs()
{
    tmp<volScalarField> tShs
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Shs",
                owner().time().timeName(),
                owner().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner().regionMesh(),
            dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& Shs = tShs();
    const scalarField& QinP = QinPrimary_.internalField();
    const scalarField& alpha = owner_.alpha().internalField();

    Shs = QinP*alpha;

    return tShs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
