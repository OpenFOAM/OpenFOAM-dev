/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionObjects::wallBoilingProperty

Description
    This function looks up the wall boiling model and collects and writes out a
    field of one of the wall boiling properties.

    Example of function object specification:
    \verbatim
    wallBoilingProperty1
    {
        type            wallBoilingProperty;
        libs            ("libmultiphaseEulerFunctionObjects.so");
        writeControl    writeTime;
        field           qEvaporative;
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallBoilingProperty | yes    |
        model    | name of the wall boiling fvModel | no   | none
        field    | name of the property       | yes        |
    \endtable

SourceFiles
    wallBoilingProperty.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoilingProperty_H
#define wallBoilingProperty_H

#include "fieldExpression.H"
#include "wallBoiling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class wallBoilingProperty Declaration
\*---------------------------------------------------------------------------*/

class wallBoilingProperty
:
    public fvMeshFunctionObject
{
private:

    //- Private Data

        //- The name of the wall boiling model
        word modelName_;

        //- The name of the property
        word fieldName_;


    // Private Member Functions

        //- Access the wall boiling model
        const fv::wallBoiling& model() const;


public:

    //- Runtime type information
    TypeName("wallBoilingProperty");


    // Constructors

        //- Construct from Time and dictionary
        wallBoilingProperty
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        wallBoilingProperty(const wallBoilingProperty&) = delete;


    //- Destructor
    virtual ~wallBoilingProperty();


    // Member Functions

        //- Read the wallBoilingProperty data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the wallBoilingProperty field
        virtual bool execute();

        //- Write the wallBoilingProperty field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wallBoilingProperty&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
