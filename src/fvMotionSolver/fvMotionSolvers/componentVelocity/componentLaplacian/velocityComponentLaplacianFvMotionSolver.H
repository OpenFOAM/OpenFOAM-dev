/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::velocityComponentLaplacianFvMotionSolver

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the given component of the motion velocity.

SourceFiles
    velocityComponentLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityComponentLaplacianFvMotionSolver_H
#define velocityComponentLaplacianFvMotionSolver_H

#include "componentVelocityMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
          Class velocityComponentLaplacianFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class velocityComponentLaplacianFvMotionSolver
:
    public componentVelocityMotionSolver,
    public fvMotionSolver
{
    // Private Data

        //- Cell-centre motion field
        mutable volScalarField cellMotionU_;

        //- Diffusivity type and data stream
        ITstream diffusivityType_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;


public:

    //- Runtime type information
    TypeName("velocityComponentLaplacian");


    // Constructors

        //- Construct from polyMesh and dictionary
        velocityComponentLaplacianFvMotionSolver
        (
            const word& name,
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        velocityComponentLaplacianFvMotionSolver
        (
            const velocityComponentLaplacianFvMotionSolver&
        ) = delete;


    //- Destructor
    ~velocityComponentLaplacianFvMotionSolver();


    // Member Functions

        //- Non-const access to the cellMotionU in order to allow changes
        //  to the boundary motion
        volScalarField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const velocityComponentLaplacianFvMotionSolver&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
