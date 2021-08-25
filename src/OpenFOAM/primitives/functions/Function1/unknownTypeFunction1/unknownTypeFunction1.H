/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::unknownTypeFunction1

Description
    Wrapper around Function1 that constructs a function for an as yet
    unknown primitive type. It stores the location of the user input and
    actually constructs the function on demand at the point at which a value is
    requested for a specific type.

SourceFiles
    unknownTypeFunction1.C
    unknownTypeFunction1Templates.C

\*---------------------------------------------------------------------------*/

#ifndef unknownTypeFunction1_H
#define unknownTypeFunction1_H

#include "Function1.H"
#include "fieldTypes.H"
#include "nil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class unknownTypeFunction1 Declaration
\*---------------------------------------------------------------------------*/

class unknownTypeFunction1
{
    // Private Classes

        //- Convenient grouping of Function1s of all primitive types
        struct functionPtr
        :
            #define PublicAutoPtrFunction1Type(Type, nullArg) \
                public autoPtr<Function1<Type>>,
            FOR_ALL_FIELD_TYPES(PublicAutoPtrFunction1Type)
            #undef InheritAutoPtrFunction1Type
            private nil
        {};


    // Private Data

        //- Name of the function
        const word name_;

        //- Top-level dictionary from which to read the function
        const dictionary& topDict_;

        //- Keyword within the top-level dictionary to the sub-dictionary in
        //  which this function is defined
        const word topDictKeyword_;

        //- Pointer to the Function1
        mutable functionPtr functionPtr_;


    // Private Member Functions

        //- Build the function of a given primitive type
        template<class Type>
        void build() const;


public:

    // Constructors

        //- Construct from an name and a dictionary
        unknownTypeFunction1(const word& name, const dictionary& dict);

        //- Disallow default bitwise copy construction
        unknownTypeFunction1(const unknownTypeFunction1&) = delete;


    //- Destructor
    ~unknownTypeFunction1();


    // Member Functions

        //- Return value as a function of scalar x
        template<class Type>
        Type value(const scalar x) const;

        //- Return value as a function of scalarField x
        template<class Type>
        tmp<Field<Type>> value(const scalarField& x) const;

        //- Integrate between two scalars
        template<class Type>
        Type integral(const scalar x1, const scalar x2) const;

        //- Integrate between two scalar fields
        template<class Type>
        tmp<Field<Type>> integral
        (
            const scalarField& x1,
            const scalarField& x2
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const unknownTypeFunction1&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "unknownTypeFunction1Templates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
