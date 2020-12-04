/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::Function2

Description
    Run-time selectable function of two variables

    with many options provided from simple constant values to complex
    functions, interpolated tabulated data etc. etc.

SourceFiles
    Function2.C
    Function2New.C

\*---------------------------------------------------------------------------*/

#ifndef Function2_H
#define Function2_H

#include "dictionary.H"
#include "tmp.H"
#include "typeInfo.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class Function2;
template<class Type> Ostream& operator<<(Ostream&, const Function2<Type>&);

/*---------------------------------------------------------------------------*\
                          Class Function2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function2
:
    public tmp<Function2<Type>>::refCount
{

protected:

    // Protected data

        //- Name of entry
        const word name_;


public:

    typedef Type returnType;

    //- Runtime type information
    TypeName("Function2")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Function2,
        dictionary,
        (
            const word& name,
            const dictionary& dict
        ),
        (name, dict)
    );


    // Constructors

        //- Construct from name
        Function2(const word& name);

        //- Copy constructor
        Function2(const Function2<Type>& f1);

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const = 0;


    //- Selector
    static autoPtr<Function2<Type>> New
    (
        const word& name,
        const dictionary& dict
    );


    //- Destructor
    virtual ~Function2();


    // Member Functions

        //- Return the name of the entry
        const word& name() const;

        //- Return value as a function of two scalars
        virtual Type value(const scalar x, const scalar y) const = 0;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const = 0;

        //- Write in dictionary format
        virtual void write(Ostream& os) const = 0;


    // Member Operators

        //- Assignment
        void operator=(const Function2<Type>&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const Function2<Type>& func
        );
};


template<class Type>
void writeEntry(Ostream& os, const Function2<Type>& f1);


/*---------------------------------------------------------------------------*\
                       Class FieldFunction2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Function2Type>
class FieldFunction2
:
    public Function2<Type>
{

public:

    // Constructors

        //- Construct from name
        FieldFunction2(const word& name);

        //- Construct and return a clone
        virtual tmp<Function2<Type>> clone() const;


    //- Destructor
    virtual ~FieldFunction2();


    // Member Functions

        //- Return value as a function of two scalars
        virtual Type value(const scalar x, const scalar y) const = 0;

        //- Return value as a function of two scalar fields
        virtual tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction2(Type)                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function2<Type>, 0);                   \
                                                                               \
    defineTemplateRunTimeSelectionTable(Function2<Type>, dictionary);


#define makeFunction2Type(SS, Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function2s::SS<Type>, 0);              \
                                                                               \
    Function2<Type>::adddictionaryConstructorToTable<Function2s::SS<Type>>     \
        addFunction2##SS##Type##ConstructorToTable_;


#define makeScalarFunction2(SS)                                                \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    Function2<scalar>::adddictionaryConstructorToTable<SS>                     \
        addFunction2##SS##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function2.C"
    #include "Constant2.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
