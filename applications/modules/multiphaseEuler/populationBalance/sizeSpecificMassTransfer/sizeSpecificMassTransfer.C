/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "sizeSpecificMassTransfer.H"
#include "groupPropertyFvScalarField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(sizeSpecificMassTransfer, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::sizeSpecificMassTransfer::readCoeffs(const dictionary& dict)
{
    if (dict.lookup<word>("populationBalance") != popBal_.name())
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change the population balance model of a " << type()
            << " model at run time" << exit(FatalIOError);
    }
}


Foam::tmp<Foam::volInternalScalarField>
Foam::fv::sizeSpecificMassTransfer::mDot() const
{
    tmp<volInternalScalarField> tmDot =
        volInternalScalarField::New
        (
            name_ + ":mDot",
            popBal_.mesh(),
            dimensionedScalar(dimensions::density/dimensions::time, scalar(0))
        );

    forAll(popBal_.fs(), groupi)
    {
        tmDot.ref() += groupMDot(groupi);
    }

    return tmDot;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::sizeSpecificMassTransfer::sizeSpecificMassTransfer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    popBal_
    (
        mesh().lookupObject<populationBalanceModel>
        (
            fvModel::coeffs(modelType, dict).lookup<word>("populationBalance")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::sizeSpecificMassTransfer::~sizeSpecificMassTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fv::sizeSpecificMassTransfer::groupi
(
    const word& fieldName
) const
{
    // At the moment we only have scalar group-fields. If we ever add
    // group-fields of other types, then this will need to cycle through them.

    if (!popBal_.mesh().foundObject<volScalarField>(fieldName)) return -1;

    const volScalarField& field =
        popBal_.mesh().lookupObject<volScalarField>(fieldName);

    if
    (
        !field.sources().table().found(name_)
     || !isA<groupPropertyFvScalarField>(field.sources().table()[name_])
    ) return -1;

    const groupPropertyFvScalarField& source =
        refCast<const groupPropertyFvScalarField>(field.sources()[name_]);

    return source.i();
}


Foam::tmp<Foam::volInternalScalarField>
Foam::fv::sizeSpecificMassTransfer::groupMDot(const label groupi) const
{
    return groupMDotByF(groupi)*popBal_.f(groupi)();
}


Foam::tmp<Foam::volInternalScalarField>
Foam::fv::sizeSpecificMassTransfer::groupS(const label groupi) const
{
    return groupMDotByF(groupi);
}


// ************************************************************************* //
