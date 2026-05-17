/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::fvMotionSolvers::velocityLaplacian

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion velocity.

SourceFiles
    velocityLaplacian_fvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityLaplacian_fvMotionSolver_H
#define velocityLaplacian_fvMotionSolver_H

#include "velocity_pointMeshMover.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

namespace fvMotionSolvers
{

/*---------------------------------------------------------------------------*\
               Class velocityLaplacian Declaration
\*---------------------------------------------------------------------------*/

class velocityLaplacian
:
    public fvMotionSolver,
    public pointMeshMovers::velocity
{
    // Private Data

        //- Cell-centre motion field
        mutable volVectorField cellMotionU_;

        //- Diffusivity type and data stream
        ITstream diffusivityType_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;


public:

    //- Runtime type information
    TypeName("velocityLaplacian");


    // Constructors

        //- Construct from fvMesh and dictionary
        velocityLaplacian(fvMesh&, const dictionary&);

        //- Construct from polyMesh and dictionary
        velocityLaplacian(const polyMesh&, const dictionary&);

        //- Disallow default bitwise copy construction
        velocityLaplacian
        (
            const velocityLaplacian&
        ) = delete;


    //- Destructor
    ~velocityLaplacian();


    // Member Functions

        //- Return reference to the cell motion velocity field
        volVectorField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return const reference to the cell motion velocity field
        const volVectorField& cellMotionU() const
        {
            return cellMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const velocityLaplacian&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMotionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
