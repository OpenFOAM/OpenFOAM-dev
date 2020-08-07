/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "UniformGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline const Foam::UniformGeometricField<Type>&
Foam::UniformGeometricField<Type>::oldTime() const
{
    return *this;
}


template<class Type>
inline const typename Foam::UniformGeometricField<Type>::Internal&
Foam::UniformGeometricField<Type>::operator()() const
{
    return *this;
}


template<class Type>
inline const typename Foam::UniformGeometricField<Type>::Internal&
Foam::UniformGeometricField<Type>::v() const
{
    return *this;
}


template<class Type>
inline typename Foam::UniformGeometricField<Type>::Internal::FieldType
Foam::UniformGeometricField<Type>::primitiveField() const
{
    return typename Internal::FieldType(this->value());
}


template<class Type>
inline typename Foam::UniformGeometricField<Type>::Boundary
Foam::UniformGeometricField<Type>::boundaryField() const
{
    return Boundary(this->value());
}


// ************************************************************************* //
