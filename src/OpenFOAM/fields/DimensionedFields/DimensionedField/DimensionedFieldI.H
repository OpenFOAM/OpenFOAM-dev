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

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const Foam::DimensionedField<Type, GeoMesh, PrimitiveField>&
Foam::DimensionedField<Type, GeoMesh, PrimitiveField>::null()
{
    return NullObjectRef<DimensionedField<Type, GeoMesh, PrimitiveField>>();
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const typename GeoMesh::Mesh&
Foam::DimensionedField<Type, GeoMesh, PrimitiveField>::mesh() const
{
    return mesh_;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const Foam::dimensionSet&
Foam::DimensionedField<Type, GeoMesh, PrimitiveField>::dimensions() const
{
    return dimensions_;
}

template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline Foam::dimensionSet&
Foam::DimensionedField<Type, GeoMesh, PrimitiveField>::dimensions()
{
    return dimensions_;
}


template<class Type, class GeoMesh, template<class> class PrimitiveField>
inline const PrimitiveField<Type>&
Foam::DimensionedField<Type, GeoMesh, PrimitiveField>::primitiveField() const
{
    return *this;
}


// ************************************************************************* //
