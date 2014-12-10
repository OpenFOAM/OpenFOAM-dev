/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "cellPointWeightWallModified.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointWeightWallModified::cellPointWeightWallModified
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label faceI
)
:
    cellPointWeight(mesh, position, cellI, faceI)
{
    // findTetrahedron or findTriangle will already have been called
    // by the cellPointWeight constructor

    if (faceI >= 0)
    {
        const polyBoundaryMesh& bm = mesh.boundaryMesh();
        label patchI = bm.whichPatch(faceI);
        if (patchI != -1)
        {
            if (isA<wallPolyPatch>(bm[patchI]))
            {
                // Apply cell centre value wall faces
                weights_[0] = 1.0;
                weights_[1] = 0.0;
                weights_[2] = 0.0;
                weights_[3] = 0.0;
            }
        }
    }
}


// ************************************************************************* //
