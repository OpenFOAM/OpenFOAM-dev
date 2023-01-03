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

#include "GeometricField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const Foam::GeometricField<Type, GeoMesh, PrimitiveField>&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::null()
{
    return NullObjectRef<GeometricField<Type, GeoMesh, PrimitiveField>>();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Internal&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::internalField() const
{
    return *this;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Internal::FieldType&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::primitiveField() const
{
    return *this;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Boundary&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::boundaryField() const
{
    return boundaryField_;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Sources&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::sources() const
{
    return sources_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Internal&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::operator()() const
{
    return *this;
}


// ************************************************************************* //
