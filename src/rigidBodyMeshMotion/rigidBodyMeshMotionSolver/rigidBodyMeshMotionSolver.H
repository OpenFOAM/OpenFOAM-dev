/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::rigidBodyMeshMotionSolver

Description
    Rigid-body mesh motion solver for fvMesh.

    Applies septernion interpolation of movement as function of distance to the
    object surface.

SourceFiles
    rigidBodyMeshMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyMeshMotionSolver_H
#define rigidBodyMeshMotionSolver_H

#include "displacementMotionSolver.H"
#include "rigidBodyMotion.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class rigidBodyMeshMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyMeshMotionSolver
:
    public motionSolver,
    public RBD::rigidBodyMotion
{
    //- Class containing the patches and point motion weighting for each body
    class bodyMesh
    {
        //- Name of the body
        const word name_;

        //- ID of the body in the RBD::rigidBodyMotion
        const label bodyID_;

        //- List of mesh patches associated with this body
        const wordReList patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;


    public:

        friend class rigidBodyMeshMotionSolver;

        bodyMesh
        (
            const polyMesh& mesh,
            const word& name,
            const label bodyID,
            const dictionary& dict
        );
    };


    // Private Data

        //- List of the bodyMeshes containing the patches and point motion
        //  weighting for each body
        PtrList<bodyMesh> bodyMeshes_;

        //- Switch for test-mode in which only the
        //  gravitational body-force is applied
        Switch test_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Ramp the forces according to the specified function and period
        autoPtr<Function1<scalar>> ramp_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        autoPtr<motionSolver> meshSolverPtr_;

        displacementMotionSolver& meshSolver_;

        //- To avoid warning from clang
        using RBD::rigidBodyMotion::write;


public:

    //- Runtime type information
    TypeName("rigidBodyMotionSolver");


    // Constructors

        //- Construct from polyMesh and dictionary
        rigidBodyMeshMotionSolver
        (
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodyMeshMotionSolver
        (
            const rigidBodyMeshMotionSolver&
        );


    //- Destructor
    ~rigidBodyMeshMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&);

        //- Write motion state information for restart
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodyMeshMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
