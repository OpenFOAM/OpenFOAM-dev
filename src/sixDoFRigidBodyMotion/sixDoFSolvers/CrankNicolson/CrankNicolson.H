/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::sixDoFSolvers::CrankNicolson

Description
    Crank-Nicolson 2nd-order time-integrator for 6DoF solid-body motion.

    The off-centering coefficients for acceleration (velocity integration) and
    velocity (position/orientation integration) may be specified but default
    values of 0.5 for each are used if they are not specified.  With the default
    off-centering this scheme is equivalent to the Newmark scheme with default
    coefficients.

    Example specification in dynamicMeshDict:
    \verbatim
    solver
    {
        type    CrankNicolson;
        aoc     0.5;    // Acceleration off-centering coefficient
        voc     0.5;    // Velocity off-centering coefficient
    }
    \endverbatim

See also
    Foam::sixDoFSolvers::Newmark

SourceFiles
    CrankNicolson.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolson_H
#define CrankNicolson_H

#include "sixDoFSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{

/*---------------------------------------------------------------------------*\
                           Class CrankNicolson Declaration
\*---------------------------------------------------------------------------*/

class CrankNicolson
:
    public sixDoFSolver
{
    // Private data

        //- Acceleration off-centering coefficient (default: 0.5)
        const scalar aoc_;

        //- Velocity off-centering coefficient (default: 0.5)
        const scalar voc_;


public:

    //- Runtime type information
    TypeName("CrankNicolson");


    // Constructors

        //- Construct from a dictionary and the body
        CrankNicolson
        (
            const dictionary& dict,
            sixDoFRigidBodyMotion& body
        );


    //- Destructor
    virtual ~CrankNicolson();


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
