/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::functionObjects::fieldsExpression

Description

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    fieldsExpression.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldsExpression_H
#define functionObjects_fieldsExpression_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class fieldsExpression Declaration
\*---------------------------------------------------------------------------*/

class fieldsExpression
:
    public fvMeshFunctionObject
{
protected:

    // Protected member data

        //- Names of fields to process
        wordList fieldNames_;

        //- Name of result fields
        word resultName_;


    // Protected member functions

        //- Set the result name
        void setResultName
        (
            const word& functionName,
            const wordList& defaultFieldNames = wordList::null()
        );

        //- Operate on the fields and store the result. Supported overload.
        template
        <
            template<class> class GeoField,
            template<class ...> class Op,
            class TypeA,
            class TypeB,
            class Enable = Op<TypeA, TypeB>
        >
        bool opAndStore
        (
            const GeoField<TypeA>& a,
            const GeoField<TypeB>& b
        );

        //- Operate on the fields and store the result. Not-supported overload.
        template
        <
            template<class> class GeoField,
            template<class ...> class Op,
            class ... Args
        >
        bool opAndStore(const Args& ...);

        //- Fold the fields expression up one place
        template
        <
            template<class> class GeoField,
            template<class ...> class Op,
            class TypeA,
            class TypeB
        >
        bool foldAB(const label i);

        //- Fold the fields expression up one place
        template
        <
            template<class> class GeoField,
            template<class ...> class Op,
            class TypeA
        >
        bool foldA(const label i);

        //- Fold the fields expression up one place
        template
        <
            template<class> class GeoField,
            template<class ...> class Op
        >
        bool fold(const label i);

        //- Calculate and store the fields expression for the given geometric
        //  field type and operation
        template<template<class> class GeoField, template<class ...> class Op>
        bool calcGeoFieldOp();

        //- Calculate and store the fields expression for the given operation
        template<template<class ...> class Op>
        bool calcOp();

        //- Perform calculation on the list of fields and return success
        virtual bool calc() = 0;


public:

    //- Runtime type information
    TypeName("fieldsExpression");


    // Constructors

        //- Construct from Time and dictionary
        fieldsExpression
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const wordList& fieldNames = wordList::null(),
            const word& resultName = word::null
        );

        //- Disallow default bitwise copy construction
        fieldsExpression(const fieldsExpression&) = delete;


    //- Destructor
    virtual ~fieldsExpression();


    // Member Functions

        //- Read the fieldsExpression data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate the result fields
        virtual bool execute();

        //- Write the result fields
        virtual bool write();

        //- Clear the result fields from the objectRegistry
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldsExpression&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldsExpressionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
