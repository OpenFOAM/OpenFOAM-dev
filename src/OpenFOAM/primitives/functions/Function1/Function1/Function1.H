/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "unitConversion.H"
#include "tmp.H"
#include "typeInfo.H"
#include "Field.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class Function1;
template<class Type> Ostream& operator<<(Ostream&, const Function1<Type>&);

/*---------------------------------------------------------------------------*\
                Class Function1s::unitConversions Declaration
\*---------------------------------------------------------------------------*/

namespace Function1s
{
struct unitConversions
{
    // Public Data

        //- Unit conversion for x-axis values
        unitConversion x;

        //- Unit conversion for result values
        unitConversion value;


    // Constructors

        //- Construct from initialiser list
        unitConversions(std::initializer_list<unitConversion>);

        //- Construct from stream
        unitConversions(Istream& is);

        //- Construct a clone
        autoPtr<unitConversions> clone() const;


    // Member Functions

        //- Update if found in the dictionary
        bool readIfPresent(const word& keyword, const dictionary&);
};
}


// Global Functions

//- Generate an error in an context where unit conversions are not supported
void assertNoConvertUnits
(
    const word& typeName,
    const Function1s::unitConversions& units,
    const dictionary& dict
);

//- Write a units entry
void writeEntry(Ostream& os, const Function1s::unitConversions& units);


// IOStream Operators

//- Read from stream
Istream& operator>>(Istream& is, Function1s::unitConversions&);

//- Write to stream
Ostream& operator<<(Ostream& os, const Function1s::unitConversions&);


/*---------------------------------------------------------------------------*\
                          Class Function1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function1
:
    public tmp<Function1<Type>>::refCount
{

protected:

    // Protected Data

        //- Name of entry
        const word name_;


public:

    typedef Type returnType;


    //- Runtime type information
    TypeName("Function1");


    // Declare runtime constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            Function1,
            dictionary,
            (
                const word& name,
                const Function1s::unitConversions& units,
                const dictionary& dict
            ),
            (name, units, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            Function1,
            Istream,
            (
                const word& name,
                const Function1s::unitConversions& units,
                Istream& is
            ),
            (name, units, is)
        );


    // Constructors

        //- Construct from name
        Function1(const word& name);

        //- Copy constructor
        Function1(const Function1<Type>& f1);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const = 0;


    // Selectors

        //- Select from dictionary
        static autoPtr<Function1<Type>> New
        (
            const word& name,
            const Function1s::unitConversions& units,
            const dictionary& dict
        );

        //- Select from dictionary
        static autoPtr<Function1<Type>> New
        (
            const word& name,
            const unitConversion& xUnits,
            const unitConversion& valueUnits,
            const dictionary& dict
        );

        //- Select from Istream
        static autoPtr<Function1<Type>> New
        (
            const word& name,
            const Function1s::unitConversions& units,
            const word& Function1Type,
            Istream& is
        );

        //- Select from Istream
        static autoPtr<Function1<Type>> New
        (
            const word& name,
            const unitConversion& xUnits,
            const unitConversion& valueUnits,
            const word& Function1Type,
            Istream& is
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

        //- Is this function guaranteed to be constant?
        virtual bool constant() const;

        //- Write data to dictionary stream
        virtual void write
        (
            Ostream& os,
            const Function1s::unitConversions&
        ) const = 0;


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

template<class Type>
void writeEntry
(
    Ostream& os,
    const Function1s::unitConversions& units,
    const Function1<Type>& f1
);

template<class Type>
void writeEntry
(
    Ostream& os,
    const unitConversion& xUnits,
    const unitConversion& valueUnits,
    const Function1<Type>& f1
);


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

#define defineFunction1(Type)                                                  \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function1<Type>, 0);                   \
    defineTemplateRunTimeSelectionTable(Function1<Type>, dictionary);          \
    defineTemplateRunTimeSelectionTable(Function1<Type>, Istream);


#define addFunction1(SS, Type)                                                 \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                          \
    typedef Function1<Type> Type##Function1;                                   \
    typedef SS<Type> Type##SS##Function1;                                      \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Type##Function1,                                                       \
        Type##SS##Function1,                                                   \
        dictionary                                                             \
    )


#define addStreamConstructableFunction1(SS, Type)                              \
                                                                               \
    addFunction1(SS, Type);                                                    \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Type##Function1,                                                       \
        Type##SS##Function1,                                                   \
        Istream                                                                \
    )


#define addNamedFunction1(SS, Type, Name)                                      \
                                                                               \
    typedef Function1<Type> Type##Function1;                                   \
    typedef SS<Type> Type##SS##Function1;                                      \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        Type##Function1,                                                       \
        Type##SS##Function1,                                                   \
        dictionary,                                                            \
        Name                                                                   \
    )


#define addScalarFunction1(SS)                                                 \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
    typedef Function1<scalar> scalarFunction1;                                 \
    addToRunTimeSelectionTable(scalarFunction1, SS, dictionary)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function1.C"
    #include "Constant.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
