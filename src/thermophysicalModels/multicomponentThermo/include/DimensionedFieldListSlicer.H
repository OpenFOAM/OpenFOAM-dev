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

Class
    Foam::DimensionedFieldListSlicer

Description
    Class to provide list slices to different parts of a dimensioned field

See also
    Foam::FieldListSlice

\*---------------------------------------------------------------------------*/

#ifndef DimensionedFieldListSlicer_H
#define DimensionedFieldListSlicer_H

#include "DimensionedField.H"
#include "FieldListSlice.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class DimensionedFieldListSlicer Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class DimensionedFieldListSlicer
{
    // Private Typedefs

        // Type of the dimensioned field
        typedef DimensionedField<Type, GeoMesh> geoFieldType;


    // Private Data

        //- Pointers to the primitive fields
        UPtrList<const Field<Type>> fields_;


public:

    // Constructors

        //- Construct null. Set later.
        inline DimensionedFieldListSlicer()
        {}

        //- Construct from a list of fields
        inline DimensionedFieldListSlicer
        (
            const UPtrList<const geoFieldType>& geoFields
        )
        {
            set(geoFields);
        }

        //- Construct from a list of fields
        inline DimensionedFieldListSlicer
        (
            const UPtrList<geoFieldType>& geoFields
        )
        :
            DimensionedFieldListSlicer
            (
                reinterpret_cast<const UPtrList<const geoFieldType>&>(geoFields)
            )
        {}


    // Member Functions

        //- Set the field pointers
        inline void set(const UPtrList<const geoFieldType>& geoFields)
        {
            fields_.resize(geoFields.size());

            forAll(geoFields, fieldi)
            {
                fields_.set
                (
                    fieldi,
                    &geoFields[fieldi].primitiveField()
                );
            }
        }

        //- Return a slice for an internal element
        inline FieldListSlice<Type> slice(const label elementi) const
        {
            return FieldListSlice<Type>(fields_, elementi);
        }
};


/*---------------------------------------------------------------------------*\
                Class DimensionedFieldListAndSlicer Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class DimensionedFieldListAndSlicer
:
    public PtrList<DimensionedField<Type, GeoMesh>>,
    private DimensionedFieldListSlicer<Type, GeoMesh>
{
    // Private Typedefs

        // Type of the dimensioned field
        typedef DimensionedField<Type, GeoMesh> geoFieldType;


public:

    // Constructors

        //- Null Constructor
        DimensionedFieldListAndSlicer()
        {}

        //- Construct with size specified
        explicit DimensionedFieldListAndSlicer(const label size)
        :
            PtrList<DimensionedField<Type, GeoMesh>>(size)
        {}

        //- Construct from arguments
        template<class ... Args>
        DimensionedFieldListAndSlicer(Args ... args)
        :
            PtrList<DimensionedField<Type, GeoMesh>>(args ...)
        {
            set();
        }


    // Member Functions

        //- Inherit the pointer list's set method
        using PtrList<DimensionedField<Type, GeoMesh>>::set;

        //- Set the field pointers in the slicer
        inline void set()
        {
            DimensionedFieldListSlicer<Type, GeoMesh>::set
            (
                reinterpret_cast<const UPtrList<const geoFieldType>&>(*this)
            );
        }

        //- Return a slice for an internal element
        using DimensionedFieldListSlicer<Type, GeoMesh>::slice;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
