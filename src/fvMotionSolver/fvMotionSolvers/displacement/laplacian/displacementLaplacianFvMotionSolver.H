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
    Foam::displacementLaplacianFvMotionSolver

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion displacement.

SourceFiles
    displacementLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementLaplacianFvMotionSolver_H
#define displacementLaplacianFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
             Class displacementLaplacianFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementLaplacianFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private Data

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Optionally read point-position field. Used only for position
        //  boundary conditions.
        mutable autoPtr<pointVectorField> pointLocation_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Frozen points (that are not on patches). -1 or points that are
        //  fixed to be at points0_ location
        label frozenPointsZone_;


public:

    //- Runtime type information
    TypeName("displacementLaplacian");


    // Constructors

        //- Construct from polyMesh and dictionary
        displacementLaplacianFvMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        displacementLaplacianFvMotionSolver
        (
            const displacementLaplacianFvMotionSolver&
        );


    //- Destructor
    ~displacementLaplacianFvMotionSolver();


    // Member Functions

        //- Return reference to the cell motion displacement field
        volVectorField& cellDisplacement()
        {
            return cellDisplacement_;
        }

        //- Return const reference to the cell motion displacement field
        const volVectorField& cellDisplacement() const
        {
            return cellDisplacement_;
        }

        //- Return reference to the diffusivity field
        motionDiffusivity& diffusivity();

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementLaplacianFvMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
