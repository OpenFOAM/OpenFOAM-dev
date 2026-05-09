/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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
    Foam::singleRigidBodyMeshMotion

Description
    Solid-body motion of the mesh specified by a run-time selectable motion
    function.  Applies SLERP interpolation of movement as function of
    distance to the object surface to move the mesh points.

SourceFiles
    singleRigidBodyMeshMotion.C

\*---------------------------------------------------------------------------*/

#ifndef singleRigidBodyMeshMotion_H
#define singleRigidBodyMeshMotion_H

#include "multiRigidBodyMeshMotion.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class singleRigidBodyMeshMotion Declaration
\*---------------------------------------------------------------------------*/

class singleRigidBodyMeshMotion
:
    public multiRigidBodyMeshMotion
{
    // Private Data

        //- The motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;


    // Private Member Functions

        //- Return the list of transformation for each body from time 0
        virtual List<septernion> transforms0() const;

        //- Move the bodies
        virtual void moveBodies();


public:

    //- Runtime type information
    TypeName("singleRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        singleRigidBodyMeshMotion
        (
            const word& name,
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        singleRigidBodyMeshMotion
        (
            const singleRigidBodyMeshMotion&
        ) = delete;


    //- Destructor
    ~singleRigidBodyMeshMotion();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singleRigidBodyMeshMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
