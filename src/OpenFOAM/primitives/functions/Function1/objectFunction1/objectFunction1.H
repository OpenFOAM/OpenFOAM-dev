/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::objectFunction1

Description
    Wrapper around Function1 that constructs a function of type associated with
    an object found in the database.

SourceFiles
    objectFunction1Templates.C

\*---------------------------------------------------------------------------*/

#ifndef objectFunction1_H
#define objectFunction1_H

#include "Function1.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class objectFunction1 Declaration
\*---------------------------------------------------------------------------*/

class objectFunction1
:
    private autoPtr<Function1<scalar>>,
    private autoPtr<Function1<vector>>,
    private autoPtr<Function1<symmTensor>>,
    private autoPtr<Function1<sphericalTensor>>,
    private autoPtr<Function1<tensor>>
{
    // Private Classes

        //- Tag structure used to deduce type in the templated constructor
        template<class PrimitiveType>
        struct type
        {};


    // Private Constructors

        //- Construct from an name and a dictionary
        template<class PrimitiveType>
        objectFunction1
        (
            const word& name,
            const dictionary& dict,
            const type<PrimitiveType>&
        );


public:

    // Selector

        template<template<class> class ObjectType>
        static autoPtr<objectFunction1> New
        (
            const word& name,
            const dictionary& dict,
            const word& objectName,
            const objectRegistry& db,
            const bool error = true
        );


    //- Destructor
    ~objectFunction1()
    {}


    // Member Functions

        //- Return value as a function of scalar x
        template<class PrimitiveType>
        PrimitiveType value(const scalar x) const;

        //- Return value as a function of scalarField x
        template<class PrimitiveType>
        tmp<Field<PrimitiveType>> value(const scalarField& x) const;

        //- Integrate between two scalars
        template<class PrimitiveType>
        PrimitiveType integral(const scalar x1, const scalar x2) const;

        //- Integrate between two scalar fields
        template<class PrimitiveType>
        tmp<Field<PrimitiveType>> integral
        (
            const scalarField& x1,
            const scalarField& x2
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const objectFunction1&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "objectFunction1Templates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
