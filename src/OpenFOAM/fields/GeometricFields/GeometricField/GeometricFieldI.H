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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
inline const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::null()
{
    return NullObjectRef<GeometricField<Type, PatchField, GeoMesh>>();
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal&
Foam::GeometricField<Type, PatchField, GeoMesh>::
internalField() const
{
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal::FieldType&
Foam::GeometricField<Type, PatchField, GeoMesh>::primitiveField() const
{
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline const typename Foam::GeometricField<Type, PatchField, GeoMesh>::
Boundary&
Foam::GeometricField<Type, PatchField, GeoMesh>::boundaryField() const
{
    return boundaryField_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline Foam::label
Foam::GeometricField<Type, PatchField, GeoMesh>::timeIndex() const
{
    return timeIndex_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
inline Foam::label&
Foam::GeometricField<Type, PatchField, GeoMesh>::timeIndex()
{
    return timeIndex_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal&
Foam::GeometricField<Type, PatchField, GeoMesh>::
operator()() const
{
    return *this;
}


// ************************************************************************* //
