/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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

#include "patchToPatchExtrapolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchToPatchExtrapolation, 0);
    defineRunTimeSelectionTable(patchToPatchExtrapolation, word);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::patchToPatchExtrapolation::update
(
    const PackedBoolList& faceCoupleds
)
{
    uncoupledFaces_ =
        selectIndices(faceCoupleds, [](const bool b) { return !b ;});

    extrapolation_ = uncoupledFaces_.size();

    reduce(extrapolation_, orOp());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatchExtrapolation::patchToPatchExtrapolation()
:
    extrapolation_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatchExtrapolation::~patchToPatchExtrapolation()
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::patchToPatchExtrapolation>
Foam::patchToPatchExtrapolation::New(const word& type)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(type);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << type << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()();
}


// ************************************************************************* //
