/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "subModelBase.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

bool Foam::subModelBase::subModelBase::inLine() const
{
    return (modelName_ != word::null);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subModelBase::subModelBase(dictionary& properties)
:
    modelName_(word::null),
    properties_(properties),
    dict_(dictionary::null),
    baseName_(word::null),
    modelType_(word::null),
    coeffDict_(dictionary::null)
{}


Foam::subModelBase::subModelBase
(
    dictionary& properties,
    const dictionary& dict,
    const word& baseName,
    const word& modelType,
    const word& dictExt
)
:
    modelName_(word::null),
    properties_(properties),
    dict_(dict),
    baseName_(baseName),
    modelType_(modelType),
    coeffDict_(dict.optionalSubDict(modelType + dictExt))
{}


Foam::subModelBase::subModelBase
(
    const word& modelName,
    dictionary& properties,
    const dictionary& dict,
    const word& baseName,
    const word& modelType
)
:
    modelName_(modelName),
    properties_(properties),
    dict_(dict),
    baseName_(baseName),
    modelType_(modelType),
    coeffDict_(dict)
{}


Foam::subModelBase::subModelBase(const subModelBase& smb)
:
    modelName_(smb.modelName_),
    properties_(smb.properties_),
    dict_(smb.dict_),
    baseName_(smb.baseName_),
    modelType_(smb.modelType_),
    coeffDict_(smb.coeffDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::subModelBase::~subModelBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::subModelBase::modelName() const
{
    return modelName_;
}


const Foam::dictionary& Foam::subModelBase::dict() const
{
    return dict_;
}


const Foam::word& Foam::subModelBase::baseName() const
{
    return baseName_;
}


const Foam::word& Foam::subModelBase::modelType() const
{
    return modelType_;
}


const Foam::dictionary& Foam::subModelBase::coeffDict() const
{
    return coeffDict_;
}


const Foam::dictionary& Foam::subModelBase::properties() const
{
    return properties_;
}


bool Foam::subModelBase::defaultCoeffs(const bool printMsg) const
{
    bool def = coeffDict_.lookupOrDefault<bool>("defaultCoeffs", false);
    if (printMsg && def)
    {
        Info<< incrIndent;
        Info<< indent << "Employing default coefficients" << endl;
        Info<< decrIndent;
    }

    return def;
}


void Foam::subModelBase::cacheFields(const bool)
{}


void Foam::subModelBase::write(Ostream& os) const
{
    os  << coeffDict_;
}


// ************************************************************************* //
