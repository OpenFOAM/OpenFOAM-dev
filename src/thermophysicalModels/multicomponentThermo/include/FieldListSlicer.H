/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::FieldListSlicer

Description
    Class to provide list slices to different parts of a field

See also
    Foam::FieldListSlice

\*---------------------------------------------------------------------------*/

#ifndef FieldListSlicer_H
#define FieldListSlicer_H

#include "Field.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class FieldListSlicer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FieldListSlicer
{
    // Private Data

        //- Slices into the field data
        List<UList<Type>> fields_;


protected:

    // Protected Member Functions

        //- Set the field pointers
        template<class Function>
        static inline void set
        (
            const Function& function,
            const label nFields,
            List<UList<Type>>& fields
        )
        {
            List<UList<Type>> tfields(nFields);
            fields.transfer(tfields);

            for (label fieldi = 0; fieldi < nFields; ++ fieldi)
            {
                fields[fieldi].shallowCopy(function(fieldi));
            }
        }

        //- Set the field pointers
        template<class Function>
        inline void set
        (
            const Function& function,
            const label nFields
        )
        {
            set(function, nFields, fields_);
        }


public:

    // Constructors

        //- Construct null. Set later.
        inline FieldListSlicer()
        {}

        //- Construct from a list of fields
        template<class Field>
        inline FieldListSlicer(const UPtrList<Field>& fields)
        {
            set(fields);
        }

        //- Construct from a sub-range of a list of fields
        template<class Field>
        inline FieldListSlicer
        (
            const UPtrList<Field>& fields,
            const label size,
            const label start
        )
        {
            set(fields, size, start);
        }


    // Member Functions

        //- Set the field pointers
        template<class Field>
        inline void set(const UPtrList<Field>& fields)
        {
            set
            (
                [&](const label fieldi) -> decltype(auto)
                {
                    return fields[fieldi];
                },
                fields.size()
            );
        }

        //- Set the field pointers to a specified sub-range
        template<class Field>
        inline void set
        (
            const UPtrList<Field>& fields,
            const label size,
            const label start
        )
        {
            set
            (
                [&](const label fieldi)
                {
                    return SubList<Type>(fields[fieldi], size, start);
                },
                fields.size()
            );
        }

        //- Return a slice for an internal element
        inline FieldListSlice<Type> slice(const label elementi) const
        {
            return FieldListSlice<Type>(fields_, elementi);
        }
};


/*---------------------------------------------------------------------------*\
                      Class FieldListAndSlicer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FieldListAndSlicer
:
    public PtrList<Field<Type>>,
    private FieldListSlicer<Type>
{
public:

    // Constructors

        //- Null Constructor
        FieldListAndSlicer()
        {}

        //- Construct with size specified
        explicit FieldListAndSlicer(const label size)
        :
            PtrList<Field<Type>>(size)
        {}

        //- Construct from arguments
        template<class ... Args>
        FieldListAndSlicer(Args ... args)
        :
            PtrList<Field<Type>>(args ...)
        {
            set();
        }


    // Member Functions

        //- Inherit the pointer list's set method
        using PtrList<Field<Type>>::set;

        //- Set the field pointers in the slicer
        inline void set()
        {
            FieldListSlicer<Type>::set(*this);
        }

        //- Return a slice for an internal element
        using FieldListSlicer<Type>::slice;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
