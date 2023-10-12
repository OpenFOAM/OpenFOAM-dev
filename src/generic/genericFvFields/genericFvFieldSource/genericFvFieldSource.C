/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "genericFvFieldSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::genericFvFieldSource<Type>::genericFvFieldSource
(
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    genericFieldBase(dict.lookup("type")),
    fvFieldSource<Type>(iF, dict),
    dict_(dict)
{}


template<class Type>
Foam::genericFvFieldSource<Type>::genericFvFieldSource
(
    const genericFvFieldSource<Type>& stf,
    const DimensionedField<Type, volMesh>& iF
)
:
    genericFieldBase(stf),
    fvFieldSource<Type>(stf, iF),
    dict_(stf.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::genericFvFieldSource<Type>::~genericFvFieldSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::genericFvFieldSource<Type>::sourceValue
(
    const fvSource& source
) const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvFieldSource"
           " (actual type " << actualTypeName() << ")"
        << "\n    on source " << source.name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic source condition."
        << abort(FatalError);

    return NullObjectRef<Field<Type>>();
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::genericFvFieldSource<Type>::internalCoeff
(
    const fvSource& source
) const
{
    FatalErrorInFunction
        << "cannot be called for a genericFvFieldSource"
           " (actual type " << actualTypeName() << ")"
        << "\n    on source " << source.name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "generic source condition."
        << abort(FatalError);

    return NullObjectRef<scalarField>();
}


template<class Type>
void Foam::genericFvFieldSource<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", actualTypeName());

    forAllConstIter(dictionary, dict_, iter)
    {
        if (iter().keyword() != "type")
        {
            iter().write(os);
        }
    }
}


// ************************************************************************* //
