/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::Function1

Description
    Run-time selectable general function of one variable

    with many options provided from simple constant values to complex
    functions, interpolated tabulated data etc. etc.

SourceFiles
    Function1.C
    Function1New.C

\*---------------------------------------------------------------------------*/

#ifndef Function1_H
#define Function1_H

#include "dictionary.H"
#include "tmp.H"
#include "typeInfo.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class Function1;
template<class Type> Ostream& operator<<(Ostream&, const Function1<Type>&);

/*---------------------------------------------------------------------------*\
                          Class Function1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function1
:
    public tmp<Function1<Type>>::refCount
{

protected:

    // Protected data

        //- Name of entry
        const word name_;


public:

    typedef Type returnType;

    //- Runtime type information
    TypeName("Function1")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Function1,
        dictionary,
        (
            const word& name,
            const dictionary& dict
        ),
        (name, dict)
    );


    // Constructors

        //- Construct from name
        Function1(const word& name);

        //- Copy constructor
        Function1(const Function1<Type>& f1);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const = 0;


    //- Selector
    static autoPtr<Function1<Type>> New
    (
        const word& name,
        const dictionary& dict
    );


    //- Destructor
    virtual ~Function1();


    // Member Functions

        //- Return the name of the entry
        const word& name() const;

        //- Return value as a function of scalar x
        virtual Type value(const scalar x) const = 0;

        //- Return value as a function of a scalar field x
        virtual tmp<Field<Type>> value(const scalarField& x) const = 0;

        //- Integrate between two scalars
        virtual Type integral(const scalar x1, const scalar x2) const = 0;

        //- Integrate between two scalar fields
        virtual tmp<Field<Type>> integral
        (
            const scalarField& x1,
            const scalarField& x2
        ) const = 0;

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const = 0;


    // Member Operators

        //- Assignment
        void operator=(const Function1<Type>&);


    // IOstream Operators

            //- Ostream Operator
            friend Ostream& operator<< <Type>
            (
                Ostream& os,
                const Function1<Type>& func
            );
};


template<class Type>
void writeEntry(Ostream& os, const Function1<Type>& f1);


/*---------------------------------------------------------------------------*\
                       Class FieldFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Function1Type>
class FieldFunction1
:
    public Function1<Type>
{

public:

    // Constructors

        //- Construct from name
        FieldFunction1(const word& name);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const;


    //- Destructor
    virtual ~FieldFunction1();


    // Member Functions

        //- Return value as a function of one scalars
        virtual Type value(const scalar x) const = 0;

        //- Return value as a function of one scalar field
        virtual tmp<Field<Type>> value(const scalarField& x) const;

        //- Integrate between two scalar values
        virtual Type integral(const scalar x1, const scalar x2) const = 0;

        //- Integrate between two scalar fields
        virtual tmp<Field<Type>> integral
        (
            const scalarField& x1,
            const scalarField& x2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction1(Type)                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function1<Type>, 0);                   \
                                                                               \
    defineTemplateRunTimeSelectionTable(Function1<Type>, dictionary);


#define makeFunction1Type(SS, Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function1s::SS<Type>, 0);              \
                                                                               \
    Function1<Type>::adddictionaryConstructorToTable<Function1s::SS<Type>>     \
    addFunction1##SS##Type##ConstructorToTable_;


#define makeNamedFunction1Type(SS, Type, Name)                                 \
                                                                               \
    Function1<Type>::adddictionaryConstructorToTable<Function1s::SS<Type>>     \
    addFunction1##Name##Type##ConstructorToTable_(#Name);


#define makeScalarFunction1(SS)                                                \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    Function1<scalar>::adddictionaryConstructorToTable<SS>                     \
        addFunction1##SS##ConstructorToTable_;


#define makeFunction1s(Type)                                                   \
                                                                               \
    template<>                                                                 \
    const char* const Foam::Tuple2<Foam::scalar, Type>::typeName               \
    (                                                                          \
        "Tuple2<scalar," #Type ">"                                             \
    );                                                                         \
                                                                               \
    makeFunction1(Type);                                                       \
    makeFunction1Type(None, Type);                                             \
    makeFunction1Type(Constant, Type);                                         \
    makeFunction1Type(Uniform, Type);                                          \
    makeFunction1Type(ZeroConstant, Type);                                     \
    makeFunction1Type(OneConstant, Type);                                      \
    makeFunction1Type(Polynomial, Type);                                       \
    makeFunction1Type(Sine, Type);                                             \
    makeFunction1Type(Square, Type);                                           \
    makeFunction1Type(Table, Type);                                            \
    makeFunction1Type(UniformTable, Type);                                     \
    makeFunction1Type(NonUniformTable, Type);                                  \
    makeNamedFunction1Type(Table, Type, tableFile);                            \
    makeFunction1Type(Scale, Type);                                            \
    makeFunction1Type(Coded, Type);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function1.C"
    #include "Constant.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
