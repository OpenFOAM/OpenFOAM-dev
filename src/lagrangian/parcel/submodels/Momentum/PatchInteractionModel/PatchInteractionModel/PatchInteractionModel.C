/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "PatchInteractionModel.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::PatchInteractionModel<CloudType>::interactionTypeNames_
(
    IStringStream
    (
        "(none rebound stick escape)"
    )()
);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::word Foam::PatchInteractionModel<CloudType>::interactionTypeToWord
(
    const interactionType& itEnum
)
{
    word it = "other";

    switch (itEnum)
    {
        case itNone:
        {
            it = "none";
            break;
        }
        case itRebound:
        {
            it = "rebound";
            break;
        }
        case itStick:
        {
            it = "stick";
            break;
        }
        case itEscape:
        {
            it = "escape";
            break;
        }
        default:
        {
        }
    }

    return it;
}


template<class CloudType>
typename Foam::PatchInteractionModel<CloudType>::interactionType
Foam::PatchInteractionModel<CloudType>::wordToInteractionType
(
    const word& itWord
)
{
    if (itWord == "none")
    {
        return itNone;
    }
    if (itWord == "rebound")
    {
        return itRebound;
    }
    else if (itWord == "stick")
    {
        return itStick;
    }
    else if (itWord == "escape")
    {
        return itEscape;
    }
    else
    {
        return itOther;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner),
    UName_("unknown_U")
{}


template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    UName_(this->coeffDict().lookupOrDefault("U", word("U")))
{}


template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    const PatchInteractionModel<CloudType>& pim
)
:
    CloudSubModelBase<CloudType>(pim),
    UName_(pim.UName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInteractionModel<CloudType>::~PatchInteractionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::word& Foam::PatchInteractionModel<CloudType>::UName() const
{
    return UName_;
}


template<class CloudType>
void Foam::PatchInteractionModel<CloudType>::info(Ostream& os)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchInteractionModelNew.C"

// ************************************************************************* //
