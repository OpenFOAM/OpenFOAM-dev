/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Collection of functions checking primitiveMesh

SourceFiles
    primitiveMeshCheck.C

\*---------------------------------------------------------------------------*/
#ifndef primitiveMeshCheck_H
#define primitiveMeshCheck_H

#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace meshCheck Declaration
\*---------------------------------------------------------------------------*/

namespace meshCheck
{
    // Topological checks

        //- Check face ordering
        bool checkUpperTriangular
        (
            const primitiveMesh& mesh,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );

        //- Check cell zip-up
        bool checkCellsZipUp
        (
            const primitiveMesh& mesh,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );

        //- Check uniqueness of face vertices
        bool checkFaceVertices
        (
            const primitiveMesh& mesh,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );

        //- Check for unused points
        bool checkPoints
        (
            const primitiveMesh& mesh,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );

        //- Check if all points on face are shared with another face.
        bool checkDuplicateFaces
        (
            const primitiveMesh& mesh,
            const label,
            const Map<label>&,
            label& nBaffleFaces,
            labelHashSet*
        );

        //- Check that shared points are in consecutive order.
        bool checkCommonOrder
        (
            const primitiveMesh& mesh,
            const label,
            const Map<label>&,
            labelHashSet*
        );

        //- Check face-face connectivity
        bool checkFaceFaces
        (
            const primitiveMesh& mesh,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );


    // Geometric checks

        //- Generate non-orthogonality field (internal faces only)
        tmp<scalarField> faceOrthogonality
        (
            const primitiveMesh& mesh,
            const vectorField& fAreas,
            const vectorField& cellCtrs
        );

        //- Generate face pyramid volume fields
        void facePyramidVolume
        (
            const primitiveMesh& mesh,
            const pointField& points,
            const vectorField& cellCtrs,
            scalarField& ownPyrVol,
            scalarField& neiPyrVol
        );

        //- Generate skewness field
        tmp<scalarField> faceSkewness
        (
            const primitiveMesh& mesh,
            const pointField& points,
            const vectorField& fCtrs,
            const vectorField& fAreas,
            const vectorField& cellCtrs
        );

        //- Generate cell openness and cell aspect ratio field
        void cellClosedness
        (
            const primitiveMesh& mesh,
            const Vector<label>& meshD,
            const vectorField& areas,
            const scalarField& vols,
            scalarField& openness,
            scalarField& aratio
        );

        //- Generate face concavity field. Returns per face the (sin of the)
        //  most concave angle between two consecutive edges
        tmp<scalarField> faceConcavity
        (
            const scalar maxSin,
            const primitiveMesh& mesh,
            const pointField& p,
            const vectorField& faceAreas
        );

        //- Generate face flatness field. Compares the individual triangles'
        //  normals against the face average normal. Between 0 (fully warped)
        //  and 1 (fully flat)
        tmp<scalarField> faceFlatness
        (
            const primitiveMesh& mesh,
            const pointField& p,
            const vectorField& fCtrs,
            const vectorField& faceAreas
        );

        //- Generate edge alignment field. Is per face the minimum aligned edge
        //  (does not use edge addressing)
        tmp<scalarField> edgeAlignment
        (
            const primitiveMesh& mesh,
            const Vector<label>& directions,
            const pointField& p
        );

        //- Generate cell determinant field
        tmp<scalarField> cellDeterminant
        (
            const primitiveMesh& mesh,
            const Vector<label>& directions,
            const vectorField& faceAreas,
            const PackedBoolList& internalOrCoupledFace
        );


        // Helpers: single face check

            //- Skewness of single face
            scalar faceSkewness
            (
                const primitiveMesh& mesh,
                const pointField& p,
                const vectorField& fCtrs,
                const vectorField& fAreas,

                const label facei,
                const point& ownCc,
                const point& neiCc
            );

            //- Skewness of single boundary face
            scalar boundaryFaceSkewness
            (
                const primitiveMesh& mesh,
                const pointField& p,
                const vectorField& fCtrs,
                const vectorField& fAreas,

                const label facei,
                const point& ownCc
            );

            //- Orthogonality of single face
            scalar faceOrthogonality
            (
                const point& ownCc,
                const point& neiCc,
                const vector& s
            );

        //- Check boundary for closedness
        bool checkClosedBoundary
        (
            const primitiveMesh& mesh,
            const scalar closedThreshold,
            const bool report = false
        );

        //- Check cells for closedness
        bool checkClosedCells
        (
            const primitiveMesh& mesh,
            const scalar closedThreshold,
            const scalar aspectThreshold,
            const bool report = false,
            labelHashSet* setPtr = nullptr,
            labelHashSet* highAspectSetPtr = nullptr,
            const Vector<label>& solutionD = Vector<label>::one
        );

        //- Check for negative face areas
        bool checkFaceAreas
        (
            const primitiveMesh& mesh,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );

        //- Check for negative cell volumes
        bool checkCellVolumes
        (
            const primitiveMesh& mesh,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );

        //- Check face pyramid volume
        bool checkFacePyramids
        (
            const primitiveMesh& mesh,
            const bool report = false,
            const scalar minPyrVol = -small,
            labelHashSet* setPtr = nullptr
        );

        //- Check face angles
        bool checkFaceAngles
        (
            const primitiveMesh& mesh,
            const bool report = false,
            const scalar maxSin = 10,    // In degrees
            labelHashSet* setPtr = nullptr
        );

        //- Check face warpage: decompose face and check ratio between
        //  magnitude of sum of triangle areas and sum of magnitude of
        //  triangle areas.
        bool checkFaceFlatness
        (
            const primitiveMesh& mesh,
            const bool report,
            const scalar warnFlatness,  // When to include in set.
            labelHashSet* setPtr
        );

        //- Check for point-point-nearness,
        //  e.g. colocated points which may be part of baffles.
        bool checkPointNearness
        (
            const primitiveMesh& mesh,
            const bool report,
            const scalar reportDistSqr,
            labelHashSet* setPtr = nullptr
        );

        //- Check edge length
        bool checkEdgeLength
        (
            const primitiveMesh& mesh,
            const bool report,
            const scalar minLenSqr,
            labelHashSet* setPtr = nullptr
        );

        //- Check for concave cells by the planes of faces
        bool checkConcaveCells
        (
            const primitiveMesh& mesh,
            const scalar planarCosAngle,
            const bool report = false,
            labelHashSet* setPtr = nullptr
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
