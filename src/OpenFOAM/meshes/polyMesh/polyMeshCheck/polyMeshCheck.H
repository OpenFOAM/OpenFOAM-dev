/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Namespace
    Foam::meshCheck

Description

SourceFiles
    polyMeshCheck.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshCheck_H
#define polyMeshCheck_H

#include "polyMesh.H"
#include "primitiveMeshCheck.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace meshCheck Declaration
\*---------------------------------------------------------------------------*/

namespace meshCheck
{
    //- Check non-orthogonality
    bool checkFaceOrthogonality
    (
        const polyMesh& mesh,
        const bool report = false,
        labelHashSet* setPtr = nullptr
    );

    //- Check face skewness
    bool checkFaceSkewness
    (
        const polyMesh& mesh,
        const bool report = false,
        labelHashSet* setPtr = nullptr
    );

    //- Check edge alignment for 1D/2D cases
    bool checkEdgeAlignment
    (
        const polyMesh& mesh,
        const bool report,
        const Vector<label>& directions,
        labelHashSet* setPtr
    );

    bool checkCellDeterminant
    (
        const polyMesh& mesh,
        const bool report,
        labelHashSet* setPtr
    );

    //- Check for face weights
    bool checkFaceWeight
    (
        const polyMesh& mesh,
        const bool report,
        const scalar minWeight = 0.05,
        labelHashSet* setPtr = nullptr
    );

    //- Check for neighbouring cell volumes
    bool checkVolRatio
    (
        const polyMesh& mesh,
        const bool report,
        const scalar minRatio = 0.01,
        labelHashSet* setPtr = nullptr
    );


    //- Check mesh topology for correctness.
    //  Returns false for no error.
    bool checkTopology(const polyMesh& mesh, const bool report = false);

    //- Check mesh geometry (& implicitly topology) for correctness.
    //  Returns false for no error.
    bool checkGeometry(const polyMesh& mesh, const bool report = false);

    //- Check mesh for correctness. Returns false for no error.
    bool checkMesh(const polyMesh& mesh, const bool report = false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
