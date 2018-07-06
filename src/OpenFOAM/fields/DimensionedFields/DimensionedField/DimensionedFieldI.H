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

template<class Type, class GeoMesh>
inline const Foam::DimensionedField<Type, GeoMesh>&
Foam::DimensionedField<Type, GeoMesh>::null()
{
    return NullObjectRef<DimensionedField<Type, GeoMesh>>();
}


template<class Type, class GeoMesh>
inline const typename GeoMesh::Mesh&
Foam::DimensionedField<Type, GeoMesh>::mesh() const
{
    return mesh_;
}


template<class Type, class GeoMesh>
inline const Foam::dimensionSet&
Foam::DimensionedField<Type, GeoMesh>::dimensions() const
{
    return dimensions_;
}

template<class Type, class GeoMesh>
inline Foam::dimensionSet&
Foam::DimensionedField<Type, GeoMesh>::dimensions()
{
    return dimensions_;
}


template<class Type, class GeoMesh>
inline const Foam::Field<Type>&
Foam::DimensionedField<Type, GeoMesh>::field() const
{
    return *this;
}

template<class Type, class GeoMesh>
inline Foam::Field<Type>&
Foam::DimensionedField<Type, GeoMesh>::field()
{
    return *this;
}


// ************************************************************************* //
