/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::dynamicMotionSolverFvMesh

Description
    The dynamicMotionSolverFvMesh

SourceFiles
    dynamicMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverFvMesh_H
#define dynamicMotionSolverFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                  Class dynamicMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverFvMesh
:
    public dynamicFvMesh
{
    // Private Data

        autoPtr<motionSolver> motionPtr_;

        //- Optional list of vectorFields to update for mesh motion
        //  For modern solvers using Uf and correctPhi to update the flux
        //  after motion it is not necessary to specify a "velocityFields" list
        velocityMotionCorrection velocityMotionCorrection_;


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicMotionSolverFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        dynamicMotionSolverFvMesh(const dynamicMotionSolverFvMesh&) = delete;


    //- Destructor
    ~dynamicMotionSolverFvMesh();


    // Member Functions

        //- Return the motionSolver
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();


        // Write

            //- Write the underlying polyMesh and other data
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write = true
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotionSolverFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
