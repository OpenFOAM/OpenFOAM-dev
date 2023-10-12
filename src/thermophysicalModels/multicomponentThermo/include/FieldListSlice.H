/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::FieldListSlice

Description
    Class to provide a list-like interface to a slice through a PtrList of
    fields

\*---------------------------------------------------------------------------*/

#ifndef FieldListSlice_H
#define FieldListSlice_H

#include "Field.H"
#include "fieldTypes.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class FieldListSlice Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class FieldListSlice
{
    // Private Data

        //- Reference to the fields
        const UPtrList<const Field<Type>>& fields_;

        //- Element index
        const label elementi_;


public:

    // Constructors

        //- Construct from a list of fields and an element index
        inline FieldListSlice
        (
            const UPtrList<const Field<Type>>& fields,
            const label elementi
        )
        :
            fields_(fields),
            elementi_(elementi)
        {}

        //- Construct from a list of fields and an element index
        inline FieldListSlice
        (
            const UPtrList<Field<Type>>& fields,
            const label elementi
        )
        :
            FieldListSlice
            (
                reinterpret_cast<const UPtrList<const Field<Type>>&>(fields),
                elementi
            )
        {}


    // Member Functions

        //- Return the size of the slice
        inline label size() const
        {
            return fields_.size();
        }


    // Member Operators

        //- Access an element
        inline const Type& operator[](const label fieldi) const
        {
            return fields_[fieldi][elementi_];
        }
};


//- Define scalarFieldListSlice, vectorFieldListSlice, etc...
#define typedefTypeFieldListSlice(Type, nullArg)                               \
    typedef FieldListSlice<Type> Type##FieldListSlice;
FOR_ALL_FIELD_TYPES(typedefTypeFieldListSlice);
#undef typedefTypeFieldListSlice


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
