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

#include "patchInteractionData.H"
#include "dictionaryEntry.H"
#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::patchInteractionData::patchInteractionData()
:
    interactionTypeName_("unknownInteractionTypeName"),
    patchName_("unknownPatch"),
    e_(0.0),
    mu_(0.0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::word& Foam::patchInteractionData::interactionTypeName() const
{
    return interactionTypeName_;
}


const Foam::word& Foam::patchInteractionData::patchName() const
{
    return patchName_;
}


Foam::scalar Foam::patchInteractionData::e() const
{
    return e_;
}


Foam::scalar Foam::patchInteractionData::mu() const
{
    return mu_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    patchInteractionData& pid
)
{
    is.check("Istream& operator>>(Istream&, patchInteractionData&)");

    const dictionaryEntry entry(dictionary::null, is);

    pid.patchName_ = entry.keyword();
    entry.lookup("type") >> pid.interactionTypeName_;
    pid.e_ = entry.lookupOrDefault<scalar>("e", 1.0);
    pid.mu_ = entry.lookupOrDefault<scalar>("mu", 0.0);

    return is;
}


// ************************************************************************* //
