/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::UniformDimensionedField

Description
    Dimensioned<Type> registered with the database as a registered IOobject
    which has the functionality of a uniform field and allows values from the
    top-level code to be passed to boundary conditions etc.

SourceFiles
    UniformDimensionedField.C

\*---------------------------------------------------------------------------*/

#ifndef UniformDimensionedField_H
#define UniformDimensionedField_H

#include "regIOobject.H"
#include "dimensionedType.H"
#include "UniformField.H"
#include "OldTimeField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class UniformDimensionedField;

template<class Type>
struct OldTimeField0Type<UniformDimensionedField<Type>>
{
    typedef UniformDimensionedField<Type> type;
};

template<class Type>
struct OldTimeOtherFieldType<UniformDimensionedField<Type>>
{
    template<template<class> class OtherPrimitiveField>
    using type = UniformDimensionedField<Type>;
};

/*---------------------------------------------------------------------------*\
                   Class UniformDimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UniformDimensionedField
:
    public regIOobject,
    public dimensioned<Type>,
    public tmp<UniformDimensionedField<Type>>::refCount,
    public OldTimeField<UniformDimensionedField<Type>>
{

public:

    // Public Typedefs

        //- Type of the non-dimensioned field
        typedef UniformField<Type> FieldType;

        //- Component type of the elements of the field
        typedef typename UniformField<Type>::cmptType cmptType;


    //- Runtime type information
    TypeName("UniformDimensionedField");


    // Constructors

        //- Construct from components. Either reads or uses supplied value.
        UniformDimensionedField(const IOobject&, const dimensioned<Type>&);

        //- Copy constructor
        UniformDimensionedField(const UniformDimensionedField<Type>&);

        //- Construct from Istream
        UniformDimensionedField(const IOobject&);


    //- Destructor
    virtual ~UniformDimensionedField();


    // Member Functions

        //- Name function provided to resolve the ambiguity between the
        //  name functions in regIOobject and dimensioned<Type>
        virtual const word& name() const
        {
            return dimensioned<Type>::name();
        }

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return a reference to the value
        Type& value();

        //- Return a const-reference to the value
        const Type& value() const;

        //- Reset the field values to the given field
        //  Equivalent to operator=
        void reset(const UniformDimensionedField<Type>&);

        //- WriteData function required for regIOobject write operation
        bool writeData(Ostream&) const;


    // Member Operators

        //- Assign to another field
        void operator==(const UniformDimensionedField<Type>&);

        //- Assign to another field
        void operator=(const UniformDimensionedField<Type>&);

        //- Assign to a dimensioned type
        void operator=(const dimensioned<Type>&);

        //- Access by index
        const Type& operator[](const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UniformDimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
