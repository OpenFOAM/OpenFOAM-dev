/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::multiValveEngineState

Description
    Writes the multi-valve engine motion state

    providing details of the position and speed of the piston and valves.

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: multiValveEngineState | yes      |
    \endtable

    Example of function object specification:
    \verbatim
    multiValveEngineState
    {
        type           multiValveEngineState;
        libs           ("libfvMeshMoversMultiValveEngine.so");
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::fvMeshMovers::multiValveEngine

SourceFiles
    multiValveEngineState.C

\*---------------------------------------------------------------------------*/

#ifndef multiValveEngineState_H
#define multiValveEngineState_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "multiValveEngine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class multiValveEngineState Declaration
\*---------------------------------------------------------------------------*/

class multiValveEngineState
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private Member Functions

        const fvMeshMovers::multiValveEngine& mve() const;


protected:

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("multiValveEngineState");


    // Constructors

        //- Construct from Time and dictionary
        multiValveEngineState
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        multiValveEngineState(const multiValveEngineState&) = delete;


    //- Destructor
    virtual ~multiValveEngineState();


    // Member Functions

        //- Read the multiValveEngineState data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the multiValveEngineState
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiValveEngineState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
