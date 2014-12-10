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
    Foam::motionSmoother

Description
    Given a displacement moves the mesh by scaling the displacement back
    until there are no more mesh errors.

    Holds displacement field (read upon construction since need boundary
    conditions) and scaling factor and optional patch number on which to
    scale back displacement.

    E.g.
    \verbatim
        // Construct iterative mesh mover.
        motionSmoother meshMover(mesh, labelList(1, patchI));

        // Set desired displacement:
        meshMover.displacement() = ..

        for (label iter = 0; iter < maxIter; iter++)
        {
            if (meshMover.scaleMesh(true))
            {
                Info<< "Successfully moved mesh" << endl;
                return true;
            }
        }
    \endverbatim

Note
    - Shared points (parallel): a processor can have points which are part of
    pp on another processor but have no pp itself (i.e. it has points
    and/or edges but no faces of pp). Hence we have to be careful when e.g.
    synchronising displacements that the value from the processor which has
    faces of pp get priority. This is currently handled in setDisplacement
    by resetting the internal displacement to zero before doing anything
    else. The combine operator used will give preference to non-zero
    values.

    - Various routines take baffles. These are sets of boundary faces that
    are treated as a single internal face. This is a hack used to apply
    movement to internal faces.

    - Mesh constraints are looked up from the supplied dictionary. (uses
    recursive lookup)

SourceFiles
    motionSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef motionSmoother_H
#define motionSmoother_H

#include "motionSmootherData.H"
#include "motionSmootherAlgo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class motionSmoother Declaration
\*---------------------------------------------------------------------------*/

class motionSmoother
:
    public motionSmootherData,
    public motionSmootherAlgo

{
    // Private class

public:

    ClassName("motionSmoother");

    // Constructors

        //- Construct from mesh, patches to work on and smoothing parameters.
        //  Reads displacement field (only boundary conditions used)
        motionSmoother
        (
            polyMesh&,
            pointMesh&,
            indirectPrimitivePatch& pp,         // 'outside' points
            const labelList& adaptPatchIDs,     // patches forming 'outside'
            const dictionary& paramDict
        );

        //- Construct from mesh, patches to work on and smoothing parameters and
        //  displacementfield (only boundary conditions used)
        motionSmoother
        (
            polyMesh&,
            indirectPrimitivePatch& pp,         // 'outside' points
            const labelList& adaptPatchIDs,     // patches forming 'outside'
            const pointVectorField&,
            const dictionary& paramDict
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
