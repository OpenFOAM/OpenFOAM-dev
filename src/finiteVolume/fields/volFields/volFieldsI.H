/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const typename
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::Internal&
Foam::GeometricField<Type, GeoMesh, PrimitiveField>::v() const
{
    static_assert(isVolMesh<GeoMesh>::value, "Only valid for volFields");
    return *this;
}


template<>
inline Foam::tmp<Foam::GeometricField<Foam::scalar, Foam::volMesh>>
Foam::GeometricField<Foam::scalar, Foam::volMesh>::component
(
    const direction
) const
{
    return *this;
}


template<>
template<template<class> class PrimitiveField>
inline void Foam::GeometricField<Foam::scalar, Foam::volMesh>::replace
(
    const direction,
    const GeometricField<scalar, volMesh, PrimitiveField>& gsf
)
{
    *this == gsf;
}


// ************************************************************************* //
