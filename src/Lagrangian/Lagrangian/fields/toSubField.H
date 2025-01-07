/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Functions to cast/convert dimensioned field references and temporaries
    based on a primitive field to equivalent types based on sub-fields

SourceFiles
    toSubField.C

\*---------------------------------------------------------------------------*/

#ifndef toSubField_H
#define toSubField_H

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a temporary sub-field from a reference to a field
template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh, SubField>> toSubField
(
    const DimensionedField<Type, GeoMesh, Field>&
);

//- Return a temporary sub-field from a reference to a sub-field. Just returns
//  the reference. No conversion necessary.
template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh, SubField>> toSubField
(
    const DimensionedField<Type, GeoMesh, SubField>&
);

//- Return a temporary sub-field from a temporary field
template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh, SubField>> toSubField
(
    const tmp<DimensionedField<Type, GeoMesh, Field>>&
);

//- Return a temporary sub-field from a temporary sub-field Just returns the
//  tmp. No conversion necessary.
template<class Type, class GeoMesh>
tmp<DimensionedField<Type, GeoMesh, SubField>> toSubField
(
    const tmp<DimensionedField<Type, GeoMesh, SubField>>&
);

//- Return a temporary sub-field from a newly constructed temporary
//  field by forwarding arguments to the sub-field's New function
template<class Type, class GeoMesh, class ... Args>
tmp<DimensionedField<Type, GeoMesh, SubField>> toSubField
(
    const Args& ... args
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "toSubFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
