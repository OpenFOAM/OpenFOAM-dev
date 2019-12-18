/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::sixDoFSolvers::Newmark

Description
    Newmark 2nd-order time-integrator for 6DoF solid-body motion.

    Reference:
    \verbatim
        Newmark, N. M. (1959).
        A method of computation for structural dynamics.
        Journal of the Engineering Mechanics Division, 85(3), 67-94.
    \endverbatim

    Example specification in dynamicMeshDict:
    \verbatim
    solver
    {
        type    Newmark;
        gamma   0.5;    // Velocity integration coefficient
        beta    0.25;   // Position integration coefficient
    }
    \endverbatim

SourceFiles
    Newmark.C

\*---------------------------------------------------------------------------*/

#ifndef Newmark_H
#define Newmark_H

#include "sixDoFSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{

/*---------------------------------------------------------------------------*\
                           Class Newmark Declaration
\*---------------------------------------------------------------------------*/

class Newmark
:
    public sixDoFSolver
{
    // Private Data

        //- Coefficient for velocity integration (default: 0.5)
        const scalar gamma_;

        //- Coefficient for position and orientation integration (default: 0.25)
        const scalar beta_;


public:

    //- Runtime type information
    TypeName("Newmark");


    // Constructors

        //- Construct from a dictionary and the body
        Newmark
        (
            const dictionary& dict,
            sixDoFRigidBodyMotion& body
        );


    //- Destructor
    virtual ~Newmark();


    // Member Functions

        //- Drag coefficient
        virtual void solve
        (
            bool firstIter,
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT,
            scalar deltaT0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sixDoFSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
