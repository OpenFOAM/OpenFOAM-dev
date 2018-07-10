/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

template<class Type>
inline Foam::label
Foam::interpolationLookUpTable<Type>::findFieldIndex
(
    const word& fieldName
) const
{
    return fieldIndices_[fieldName];
}


template<class Type>
inline const Foam::List<Foam::dictionary>&
Foam::interpolationLookUpTable<Type>::output() const
{
    return output_;
}


template<class Type>
inline const Foam::List<Foam::dictionary>&
Foam::interpolationLookUpTable<Type>::entries() const
{
    return entries_;
}


template<class Type>
inline const Foam::List<Foam::scalar>&
Foam::interpolationLookUpTable<Type>::min() const
{
    return min_;
}


template<class Type>
inline const Foam::List<Foam::label>&
Foam::interpolationLookUpTable<Type>::dim() const
{
    return dim_;
}


template<class Type>
inline const Foam::List<Foam::scalar>&
Foam::interpolationLookUpTable<Type>::delta() const
{
    return delta_;
}


template<class Type>
inline const Foam::List<Foam::scalar>&
Foam::interpolationLookUpTable<Type>::max() const
{
    return max_;
}


template<class Type>
inline Foam::word Foam::interpolationLookUpTable<Type>::tableName() const
{
    return fileName_.name();
}


// ************************************************************************* //
