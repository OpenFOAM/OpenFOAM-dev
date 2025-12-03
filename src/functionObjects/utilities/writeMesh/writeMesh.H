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
    Foam::functionObjects::writeMesh

Description
    Writes changing mesh data (e.g., points, faces, zones) at times other than
    the write times defined in the controlDict

    Example of function object specification:
    \verbatim
    writeMesh1
    {
        type            writeMesh;
        libs            ("libutilityFunctionObjects.so");
        writeControl    adjustableRunTime;
        writeInterval   0.002;
    }
    \endverbatim

Usage
    \table
        Property | Description             | Required | Default value
        type     | type name: writeMesh    | yes      |
        region   | the name of the region  | no       | default region
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    writeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef writeMesh_functionObject_H
#define writeMesh_functionObject_H

#include "regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class writeMesh Declaration
\*---------------------------------------------------------------------------*/

class writeMesh
:
    public regionFunctionObject
{
public:

    //- Runtime type information
    TypeName("writeMesh");


    // Constructors

        //- Inherit constructors
        using regionFunctionObject::regionFunctionObject;


    //- Destructor
    virtual ~writeMesh();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Return true if the functionObject should be executed at the start
        virtual bool executeAtStart() const;

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
