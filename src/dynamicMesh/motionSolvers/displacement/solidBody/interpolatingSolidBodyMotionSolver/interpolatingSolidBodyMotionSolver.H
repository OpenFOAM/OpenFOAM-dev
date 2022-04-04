/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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
    Foam::interpolatingSolidBodyMotionSolver

Description
    Solid-body motion of the mesh specified by a run-time selectable motion
    function.  Applies SLERP interpolation of movement as function of
    distance to the object surface to move the mesh points.

SourceFiles
    interpolatingSolidBodyMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatingSolidBodyMotionSolver_H
#define interpolatingSolidBodyMotionSolver_H

#include "points0MotionSolver.H"
#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class interpolatingSolidBodyMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class interpolatingSolidBodyMotionSolver
:
    public points0MotionSolver
{
    // Private Data

        //- The motion control function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        wordReList patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Center of gravity read from dictionary
        vector CofG_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Current interpolation scale (1 at patches, 0 at distance_)
        pointScalarField scale_;


public:

    //- Runtime type information
    TypeName("interpolatingSolidBody");


    // Constructors

        //- Construct from polyMesh and dictionary
        interpolatingSolidBodyMotionSolver
        (
            const word& name,
            const polyMesh&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        interpolatingSolidBodyMotionSolver
        (
            const interpolatingSolidBodyMotionSolver&
        );


    //- Destructor
    ~interpolatingSolidBodyMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interpolatingSolidBodyMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
