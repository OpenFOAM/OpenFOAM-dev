/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::forwardOrAssignFieldMapper

Description
    Forward field mapper that permits assigning a value on unmapped elements

\*---------------------------------------------------------------------------*/

#ifndef forwardOrAssignFieldMapper_H
#define forwardOrAssignFieldMapper_H

#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class forwardOrAssignFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class forwardOrAssignFieldMapper
:
    public fieldMapper
{
    // Private Data

        //- Addressing from new back to old
        const labelUList& addressing_;

        //- Are any values unmapped?
        bool hasUnmapped_;


    // Private Member Functions

        //- Generate an error when parts of the field are unmapped, but no
        //  value has been provided for them
        virtual void unmappedError() const;

        //- Map a field in-place
        template<class Type>
        void map(Field<Type>& f, const Field<Type>& mapF) const;

        //- Map a field and return the result
        template<class Type>
        tmp<Field<Type>> map(const Field<Type>& f) const;

        //- Map or assign a field in-place
        template<class Type>
        void mapOrAssign
        (
            Field<Type>& f,
            const Field<Type>& mapF,
            const Type& unmappedVal
        ) const;

        //- Map or assign a field and return the result
        template<class Type>
        tmp<Field<Type>> mapOrAssign
        (
            const Field<Type>& mapF,
            const Type& unmappedVal
        ) const;

        //- Map or assign a field in-place
        template<class Type>
        void mapOrAssign
        (
            Field<Type>& f,
            const Field<Type>& mapF,
            const FieldFunctor<Type>& unmappedFunc
        ) const;

        //- Map or assign a field and return the result
        template<class Type>
        tmp<Field<Type>> mapOrAssign
        (
            const Field<Type>& mapF,
            const FieldFunctor<Type>& unmappedFunc
        ) const;


public:

    // Constructors

        //- Construct given addressing
        forwardOrAssignFieldMapper(const labelUList& addressing)
        :
            addressing_(addressing),
            hasUnmapped_(!addressing.empty() && min(addressing) < 0)
        {}


    //- Destructor
    virtual ~forwardOrAssignFieldMapper()
    {}


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_MAP_OPERATOR, );

        //- Map a label field
        DEFINE_FIELD_MAPPER_MAP_OPERATOR(label, );

        //- Map or assign a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR,);

        //- Map or assign a label field
        DEFINE_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR(label,);
};


/*---------------------------------------------------------------------------*\
              Class forwardOrAssignPatchFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class forwardOrAssignPatchFieldMapper
:
    public forwardOrAssignFieldMapper
{
    // Private Data

        //- Field base type name
        const word& fieldBaseTypeName_;

        //- Field type name
        const word& fieldTypeName_;

        //- Patch name
        const word& patchName_;

        //- Internal patch name
        const word& internalFieldName_;


    // Private Member Functions

        //- Generate an error when parts of the field are unmapped, but no
        //  value has been provided for them
        virtual void unmappedError() const;


public:

    // Constructors

        //- Construct given a field and addressing
        template<class FieldType>
        forwardOrAssignPatchFieldMapper
        (
            const FieldType& field,
            const labelUList& addressing
        )
        :
            forwardOrAssignFieldMapper(addressing),
            fieldBaseTypeName_(FieldType::typeName),
            fieldTypeName_(field.type()),
            patchName_(field.patch().name()),
            internalFieldName_(field.internalField().name())
        {}


    //- Destructor
    virtual ~forwardOrAssignPatchFieldMapper()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
