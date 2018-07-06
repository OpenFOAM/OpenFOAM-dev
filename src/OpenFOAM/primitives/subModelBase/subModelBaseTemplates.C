/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::subModelBase::getBaseProperty
(
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;

    if (properties_.found(baseName_))
    {
        const dictionary& baseDict = properties_.subDict(baseName_);
        baseDict.readIfPresent(entryName, result);
    }

    return result;
}


template<class Type>
void Foam::subModelBase::getBaseProperty
(
    const word& entryName,
    Type& value
) const
{
    if (properties_.found(baseName_))
    {
        const dictionary& baseDict = properties_.subDict(baseName_);
        baseDict.readIfPresent(entryName, value);
    }
}


template<class Type>
void Foam::subModelBase::setBaseProperty
(
    const word& entryName,
    const Type& value
)
{
    if (properties_.found(baseName_))
    {
        dictionary& baseDict = properties_.subDict(baseName_);
        baseDict.add(entryName, value, true);
    }
    else
    {
        properties_.add(baseName_, dictionary());
        properties_.subDict(baseName_).add(entryName, value);
    }
}


template<class Type>
void Foam::subModelBase::getModelProperty
(
    const word& entryName,
    Type& value
) const
{
    if (properties_.found(baseName_))
    {
        const dictionary& baseDict = properties_.subDict(baseName_);

        if (inLine() && baseDict.found(modelName_))
        {
            baseDict.subDict(modelName_).readIfPresent(entryName, value);
        }
        else if (baseDict.found(modelType_))
        {
            baseDict.subDict(modelType_).readIfPresent(entryName, value);
        }
    }
}


template<class Type>
Type Foam::subModelBase::getModelProperty
(
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;
    getModelProperty(entryName, result);
    return result;
}


template<class Type>
void Foam::subModelBase::setModelProperty
(
    const word& entryName,
    const Type& value
)
{
    if (properties_.found(baseName_))
    {
        dictionary& baseDict = properties_.subDict(baseName_);

        if (inLine())
        {
            if (baseDict.found(modelName_))
            {
                baseDict.subDict(modelName_).add(entryName, value, true);
            }
            else
            {
                baseDict.add(modelName_, dictionary());
                baseDict.subDict(modelName_).add(entryName, value, true);
            }
        }
        else
        {
            if (baseDict.found(modelType_))
            {
                baseDict.subDict(modelType_).add(entryName, value, true);
            }
            else
            {
                baseDict.add(modelType_, dictionary());
                baseDict.subDict(modelType_).add(entryName, value, true);
            }
        }
    }
    else
    {
        properties_.add(baseName_, dictionary());

        if (inLine())
        {
            properties_.subDict(baseName_).add(modelName_, dictionary());
            properties_.subDict(baseName_).subDict(modelName_).add
            (
                entryName,
                value
            );
        }
        else
        {
            properties_.subDict(baseName_).add(modelType_, dictionary());
            properties_.subDict(baseName_).subDict(modelType_).add
            (
                entryName,
                value
            );
        }
    }
}


// ************************************************************************* //
