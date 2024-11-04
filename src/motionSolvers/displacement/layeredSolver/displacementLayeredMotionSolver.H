/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::displacementLayeredMotionSolver

Description
    Interpolating motion solver for extruded/layered meshes

    Uses a topological walk from points along edges to points from one patch to
    the opposite, accumulating distance and propagating displacement and then
    again in the return direction and uses the distances from each point to the
    two opposite patches to linearly interpolate the displacements of
    these two patches to the each point.

    This motion solver is very simple to configure, only the two opposite
    patches between which the mesh is layered need be specified as a pair of
    patch names following the \c oppositePatches keyword, e.g.

    \verbatim
        mover
        {
            type            motionSolver;

            libs            ("libfvMotionSolvers.so");

            motionSolver    displacementLayered;

            oppositePatches (top bottom);
        }
    \endverbatim

SourceFiles
    displacementLayeredMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementLayeredMotionSolver_H
#define displacementLayeredMotionSolver_H

#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class displacementLayeredMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementLayeredMotionSolver
:
    public displacementMotionSolver
{
    // Private Member Data

        //- Names of the opposite patches of the layered mesh
        Pair<word> oppositePatchNames_;

        //- Indices of the opposite patches of the layered mesh
        Pair<label> oppositePatches_;


    // Private Member Functions

        //- Walk across the layers from the startPatch to the opposite patch
        //  accumulating the distance and setting the displacement
        void walkLayers
        (
            const polyPatch& startPatch,
            scalarField& distance,
            vectorField& displacement
        ) const;


public:

    //- Runtime type information
    TypeName("displacementLayered");


    // Constructors

        //- Construct from polyMesh and dictionary
        displacementLayeredMotionSolver
        (
            const word& name,
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        displacementLayeredMotionSolver
        (
            const displacementLayeredMotionSolver&
        ) = delete;


    //- Destructor
    ~displacementLayeredMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology (not implemented)
        virtual void topoChange(const polyTopoChangeMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const displacementLayeredMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
