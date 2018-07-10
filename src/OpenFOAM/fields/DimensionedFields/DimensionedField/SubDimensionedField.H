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

Class
    Foam::SubDimensionedField

Description
    SubDimensionedField is a DimensionedField obtained as a section of another
    DimensionedField.

    Thus it is itself unallocated so that no storage is allocated or
    deallocated during its use.  To achieve this behaviour,
    SubDimensionedField is derived from SubField rather than Field.

SourceFiles
    SubDimensionedFieldI.H

\*---------------------------------------------------------------------------*/

#ifndef SubDimensionedField_H
#define SubDimensionedField_H

#include "Field.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SubDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class SubDimensionedField
:
    public regIOobject,
    public SubField<Type>
{

public:

    // Public typedefs

        typedef typename GeoMesh::Mesh Mesh;
        typedef typename Field<Type>::cmptType cmptType;


    // Constructors

        //- Construct from a SubField
        inline SubDimensionedField
        (
            const SubField<Type>& slist
        );

        //- Construct from a UList and size
        inline SubDimensionedField
        (
            const UList<Type>& list,
            const label subSize
        );

        //- Construct from a UList start and end indices
        inline SubDimensionedField
        (
            const UList<Type>& list,
            const label subSize,
            const label startIndex
        );

        //- Construct as copy
        inline SubDimensionedField
        (
            const SubDimensionedField<cmptType, GeoMesh>& sfield
        );


    // Member functions

        //- Return a null SubDimensionedField
        static inline const SubDimensionedField<Type, GeoMesh>& null();

        //- Return a component field of the field
        inline tmp<DimensionedField<cmptType, GeoMesh>> component
        (
            const direction
        ) const;

        //- Return the field transpose (only defined for second rank tensors)
        tmp<DimensionedField<Type, GeoMesh>> T() const;


    // Member operators

        //- Assignment
        inline void operator=(const SubDimensionedField<Type, GeoMesh>&);

        //- Allow cast to a const DimensionedField<Type, GeoMesh>&
        inline operator const DimensionedField<Type, GeoMesh>&() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubDimensionedFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
