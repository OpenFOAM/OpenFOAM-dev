/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::sixDoFRigidBodyMotionSolver

Description
    6-DoF solid-body mesh motion solver for an fvMesh.

    Applies SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    sixDoFRigidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionSolver_H
#define sixDoFRigidBodyMotionSolver_H

#include "displacementMotionSolver.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sixDoFRigidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionSolver
:
    public displacementMotionSolver,
    public sixDoFRigidBodyMotion
{
    // Private Data

        wordReList patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

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

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and dictionary
        sixDoFRigidBodyMotionSolver
        (
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sixDoFRigidBodyMotionSolver
        (
            const sixDoFRigidBodyMotionSolver&
        );


    //- Destructor
    ~sixDoFRigidBodyMotionSolver();


    // Member Functions

        //- Return the six DoF motion object
        const sixDoFRigidBodyMotion& motion() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Write motion state information for restart
        virtual bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sixDoFRigidBodyMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
