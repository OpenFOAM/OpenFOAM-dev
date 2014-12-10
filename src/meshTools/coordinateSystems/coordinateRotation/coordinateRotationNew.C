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

#include "coordinateRotation.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateRotation> Foam::coordinateRotation::New
(
    const dictionary& dict, const objectRegistry& obr
)
{
    if (debug)
    {
        Pout<< "coordinateRotation::New(const dictionary&) : "
            << "constructing coordinateRotation"
            << endl;
    }

    word rotType = dict.lookup("type");

    objectRegistryConstructorTable::iterator cstrIter =
        objectRegistryConstructorTablePtr_->find(rotType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "coordinateRotation::New"
            "("
            "   const dictionary&, "
            "   const objectRegistry& "
            ")",
            dict
        )   << "Unknown coordinateRotation type "
            << rotType << nl << nl
            << "Valid coordinateRotation types are :" <<  nl
            << "[default: axes ]"
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateRotation>(cstrIter()(dict, obr));
}


Foam::autoPtr<Foam::coordinateRotation> Foam::coordinateRotation::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        Pout<< "coordinateRotation::New(const dictionary&) : "
            << "constructing coordinateRotation"
            << endl;
    }

    word rotType = dict.lookup("type");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(rotType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "coordinateRotation::New"
            "("
            "   const dictionary&, "
            ")",
            dict
        )   << "Unknown coordinateRotation type "
            << rotType << nl << nl
            << "Valid coordinateRotation types are :" <<  nl
            << "[default: axes ]"
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateRotation>(cstrIter()(dict));
}

// ************************************************************************* //
