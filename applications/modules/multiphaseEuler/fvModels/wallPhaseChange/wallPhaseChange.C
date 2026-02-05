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

#include "wallPhaseChange.H"
#include "zeroFixedValueFvPatchFields.H"
#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(wallPhaseChange, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::fv::wallPhaseChange::mDotBoundaryTypes
(
    const word& activeType
) const
{
    wordList boundaryTypes
    (
        mesh().boundary().size(),
        zeroFixedValueFvPatchScalarField::typeName
    );

    forAll(fluid_.mesh().boundary(), patchi)
    {
        const bool is1Active =
            isA<alphatPhaseChangeWallFunctionFvPatchScalarField>
            (
                alphats_.first().boundaryField()[patchi]
            );
        const bool is2Active =
            isA<alphatPhaseChangeWallFunctionFvPatchScalarField>
            (
                alphats_.second().boundaryField()[patchi]
            );

        if (is1Active != is2Active)
        {
            FatalErrorInFunction
                << "The field "
                << (is1Active ? alphats_.first() : alphats_.second()).name()
                << " has a phase change wall function on patch "
                << mesh().boundary()[patchi].name() << " but "
                << (is2Active ? alphats_.first() : alphats_.second()).name()
                << " does not" << exit(FatalError);
        }

        if (is1Active)
        {
            boundaryTypes[patchi] = activeType;
        }
    }

    return boundaryTypes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::wallPhaseChange::wallPhaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const wordList& species
)
:
    phaseChange(name, modelType, mesh, dict, species),
    fluid_(mesh().lookupObject<phaseSystem>(phaseSystem::propertiesName)),
    phases_
    (
        fluid_.phases()[phaseNames().first()],
        fluid_.phases()[phaseNames().second()]
    ),
    alphats_
    (
        mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alphat", phases_.first().name())
        ),
        mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alphat", phases_.second().name())
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::wallPhaseChange::isPatchActive(const label patchi) const
{
    return
        isA<alphatPhaseChangeWallFunctionFvPatchScalarField>
        (
            alphats_.first().boundaryField()[patchi]
        )
     || isA<alphatPhaseChangeWallFunctionFvPatchScalarField>
        (
            alphats_.second().boundaryField()[patchi]
        );
}


// ************************************************************************* //
