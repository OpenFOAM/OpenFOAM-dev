/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::functionObjects::meshToMeshAdjustTimeStepFunctionObject

Description
    Adjusts time-step for meshToMesh mapping.

SourceFiles
    meshToMeshAdjustTimeStepFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_meshToMeshAdjustTimeStepFunctionObject_H
#define functionObjects_meshToMeshAdjustTimeStepFunctionObject_H

#include "fvMeshFunctionObject.H"
#include "fvMeshTopoChangersMeshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class meshToMeshAdjustTimeStepFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class meshToMeshAdjustTimeStepFunctionObject
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Reference to the meshToMesh fvMeshTopoChanger
        const fvMeshTopoChangers::meshToMesh& meshToMesh_;


public:

    //- Runtime type information
    TypeName("meshToMeshAdjustTimeStep");


    // Constructors

        //- Construct from components
        meshToMeshAdjustTimeStepFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        meshToMeshAdjustTimeStepFunctionObject
        (
            const meshToMeshAdjustTimeStepFunctionObject&
        ) = delete;


    // Destructor
    virtual ~meshToMeshAdjustTimeStepFunctionObject();


    // Member Functions

        //-  Read and reset the timeStep Function1
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Return the time to the next write
        virtual scalar timeToNextWrite();

        //- Do nothing
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshToMeshAdjustTimeStepFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
