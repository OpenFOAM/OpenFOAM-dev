/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::fieldMapper

Description
    Abstract base class for field mapping

\*---------------------------------------------------------------------------*/

#ifndef fieldMapper_H
#define fieldMapper_H

#include "Field.H"
#include "fieldTypes.H"
#include "fieldMapperM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fieldMapper Declaration
\*---------------------------------------------------------------------------*/

class fieldMapper
{
public:

    // Public Classes

        //- Class used to lazily evaluate fields
        template<class Type>
        class FieldFunctor;

        //- Class used to lazily evaluate field-generating operators
        template<class Type, class FieldOp>
        class FieldOpFunctor;


private:

    // Private Typedefs

        //- Alias for tmp<Field<Type>> if the other argument passed is an
        //  operator (i.e., not a field). Disambiguates calls that take
        //  operators from those that take field arguments.
        template<class Type, class FieldOp>
        using TmpFieldTypeIfFieldOp =
            typename std::enable_if
            <
                !std::is_base_of<Field<Type>, FieldOp>::value,
                tmp<Field<Type>>
            >::type;


    // Private Member Functions

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

        //- Null constructor
        fieldMapper()
        {}


    //- Destructor
    virtual ~fieldMapper()
    {}


    // Member Operators

        //- Map a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_MAP_OPERATOR, = 0);

        //- Map a label field
        DEFINE_FIELD_MAPPER_MAP_OPERATOR(label, = 0);

        //- Map a temporary field
        template<class Type>
        void operator()(Field<Type>& f, const tmp<Field<Type>>& tmapF) const;

        //- Map a temporary field
        template<class Type>
        tmp<Field<Type>> operator()(const tmp<Field<Type>>& tmapF) const;

        //- Map or assign a field
        FOR_ALL_FIELD_TYPES(DEFINE_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR,);

        //- Map or assign a label field
        DEFINE_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR(label,);

        //- Map or assign a field from an operator in-place
        template<class Type, class FieldOp>
        void operator()
        (
            Field<Type>& f,
            const Field<Type>& mapF,
            const FieldOp& unmappedOp
        ) const;

        //- Map or assign a field from an operator and return the result
        template<class Type, class FieldOp>
        TmpFieldTypeIfFieldOp<Type, FieldOp> operator()
        (
            const Field<Type>& mapF,
            const FieldOp& unmappedOp
        ) const;
};


/*---------------------------------------------------------------------------*\
                  Class fieldMapper::FieldFunctor Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fieldMapper::FieldFunctor
{
public:

    // Constructors

        //- Construct null
        FieldFunctor()
        {}


    //- Destructor
    virtual ~FieldFunctor()
    {}


    // Member Operators

        //- Evaluate the field
        virtual tmp<Field<Type>> operator()() const = 0;
};


/*---------------------------------------------------------------------------*\
                  Class fieldMapper::FieldOpFunctor Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class FieldOp>
class fieldMapper::FieldOpFunctor
:
    public FieldFunctor<Type>
{
    // Private Member Data

        //- The operator
        FieldOp op_;


public:

    // Constructors

        //- Construct from an operator
        FieldOpFunctor(const FieldOp& op)
        :
            FieldFunctor<Type>(),
            op_(op)
        {}


    //- Destructor
    virtual ~FieldOpFunctor()
    {}


    // Member Operators

        //- Evaluate the field
        virtual tmp<Field<Type>> operator()() const
        {
            return op_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
