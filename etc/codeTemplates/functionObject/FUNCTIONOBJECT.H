/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    Foam::functionObjects::FUNCTIONOBJECT

Description
    This function object...

    Example of function object specification:
    \verbatim
    FUNCTIONOBJECT1
    {
        type           FUNCTIONOBJECT;
        libs ("libFUNCTIONOBJECTFunctionObject.so");
        ...
        wordData       someWord;
        scalarData     1.0;
        labelData      1;
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: FUNCTIONOBJECT | yes      |
        wordData     | some word option...       | no       | defaultWord
        scalarData   | some scalar value...      | yes      |
        labelData    | some label value...       | yes      |
    \endtable

SourceFiles
    FUNCTIONOBJECT.C

\*---------------------------------------------------------------------------*/

#ifndef FUNCTIONOBJECT_H
#define FUNCTIONOBJECT_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class FUNCTIONOBJECT Declaration
\*---------------------------------------------------------------------------*/

class FUNCTIONOBJECT
:
    public fvMeshFunctionObject
{
    // Private data

        //- word
        word wordData_;

        //- scalar
        scalar scalarData_;

        //- label
        label labelData_;


    // Private Member Functions


        //- Disallow default bitwise copy construct
        FUNCTIONOBJECT(const FUNCTIONOBJECT&);

        //- Disallow default bitwise assignment
        void operator=(const FUNCTIONOBJECT&);


public:

    //- Runtime type information
    TypeName("FUNCTIONOBJECT");


    // Constructors

        //- Construct from Time and dictionary
        FUNCTIONOBJECT
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~FUNCTIONOBJECT();


    // Member Functions

        //- Read the FUNCTIONOBJECT data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the FUNCTIONOBJECT
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
