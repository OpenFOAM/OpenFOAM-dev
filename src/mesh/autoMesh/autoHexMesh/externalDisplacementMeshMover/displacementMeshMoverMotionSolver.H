/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::displacementMeshMoverMotionSolver

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion displacement.

SourceFiles
    displacementMeshMoverMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementMeshMoverMotionSolver_H
#define displacementMeshMoverMotionSolver_H

#include "displacementMotionSolver.H"
#include "externalDisplacementMeshMover.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class displacementMeshMoverMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementMeshMoverMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        mutable autoPtr<externalDisplacementMeshMover> meshMoverPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        displacementMeshMoverMotionSolver
        (
            const displacementMeshMoverMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const displacementMeshMoverMotionSolver&);


public:

    //- Runtime type information
    TypeName("displacementMeshMover");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        displacementMeshMoverMotionSolver(const polyMesh&, const IOdictionary&);


    //- Destructor
    ~displacementMeshMoverMotionSolver();


    // Member Functions

        externalDisplacementMeshMover& meshMover() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
