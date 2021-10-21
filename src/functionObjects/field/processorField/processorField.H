/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::functionObjects::processorField

Description
    Writes a scalar field whose value is the local processor ID.  The output
    field name is 'processorID'.

    Example of function object specification:
    \verbatim
    processorField1
    {
        type        processorField;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: processorField | yes       |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    processorField.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_processorField_H
#define functionObjects_processorField_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class processorField Declaration
\*---------------------------------------------------------------------------*/

class processorField
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("processorField");


    // Constructors

        //- Construct from Time and dictionary
        processorField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        processorField(const processorField&) = delete;


    //- Destructor
    virtual ~processorField();


    // Member Functions

        //- Read the input data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the processorID field
        virtual bool execute();

        //- Write the processorID field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const processorField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
