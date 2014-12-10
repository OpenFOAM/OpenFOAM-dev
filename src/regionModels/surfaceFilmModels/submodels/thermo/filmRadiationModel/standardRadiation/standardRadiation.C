/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    standardRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardRadiation::standardRadiation
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
    ),
    QrNet_
    (
        IOobject
        (
            "QrNet",
            owner.time().timeName(),
            owner.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        owner.regionMesh(),
        dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    beta_(readScalar(coeffDict_.lookup("beta"))),
    kappaBar_(readScalar(coeffDict_.lookup("kappaBar")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardRadiation::~standardRadiation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardRadiation::correct()
{
    // Transfer Qr from primary region
    QinPrimary_.correctBoundaryConditions();
}


tmp<volScalarField> standardRadiation::Shs()
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
    const scalarField& delta = owner_.delta().internalField();
    const scalarField& alpha = owner_.alpha().internalField();

    Shs = beta_*QinP*alpha*(1.0 - exp(-kappaBar_*delta));

    // Update net Qr on local region
    QrNet_.internalField() = QinP - Shs;
    QrNet_.correctBoundaryConditions();

    return tShs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
