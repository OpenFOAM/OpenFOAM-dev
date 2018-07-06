/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

        void setResultName
        (
            const word& typeName,
            const wordList& defaultArg = wordList::null()
        );

        //- Call 'calcFieldType' for the given functionObject
        //  for 'volField' and 'surfaceField' field types
        template<class Type, class FOType>
        bool calcFieldTypes(FOType& fo);

        //- Call 'calcFieldTypes' for the given 'Type' and functionObject
        template<class Type, class FOType>
        bool calcType(FOType& fo);

        //- Call 'calcType' for the given functionObject
        //  for each primitive type
        template<class FOType>
        bool calcAllTypes(FOType& fo);

        virtual bool calc() = 0;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fieldsExpression(const fieldsExpression&);

        //- Disallow default bitwise assignment
        void operator=(const fieldsExpression&);


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


    //- Destructor
    virtual ~fieldsExpression();


    // Member Functions

        //- Read the fieldsExpression data
        virtual bool read(const dictionary&);

        //- Calculate the result fields
        virtual bool execute();

        //- Write the result fields
        virtual bool write();

        //- Clear the result fields from the objectRegistry
        virtual bool clear();
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
