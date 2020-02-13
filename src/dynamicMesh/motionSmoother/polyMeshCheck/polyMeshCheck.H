/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Description
    PolyMesh checking routines. Checks various criteria for a mesh and supplied
    geometry, with the mesh only used for topology. Coupled patch aware (i.e.,
    coupled faces are treated as internal).

SourceFiles
    polyMeshCheck.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshCheck_H
#define polyMeshCheck_H

#include "pointFields.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace polyMeshCheck
{

//- Check face non-orthogonality
bool checkFaceDotProduct
(
    const bool report,
    const scalar orthWarn,
    const polyMesh&,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
);

//- Check face pyramid volumes
bool checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    const polyMesh&,
    const vectorField& cellCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet*
);

//- Check face tetrahedra volumes
bool checkFaceTets
(
    const bool report,
    const scalar minPyrVol,
    const polyMesh&,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet*
);

//- Check face skewness
bool checkFaceSkewness
(
    const bool report,
    const scalar internalSkew,
    const scalar boundarySkew,
    const polyMesh& mesh,
    const pointField& points,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
);

//- Check interpolation weights (0.5 for regular mesh)
bool checkFaceWeights
(
    const bool report,
    const scalar warnWeight,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
);

//- Cell volume ratio of neighbouring cells (1 for regular mesh)
bool checkVolRatio
(
    const bool report,
    const scalar warnRatio,
    const polyMesh& mesh,
    const scalarField& cellVolumes,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
);

//- Check convexity of angles in a face. See primitiveMesh for explanation.
bool checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
);

// Check the difference between normals of individual face-triangles (from
// face-centre decomposition) and the cell-cell centre vector
bool checkFaceTwist
(
    const bool report,
    const scalar minTwist,
    const polyMesh&,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
);

//- Check consecutive face-triangle (from face-centre decomposition) normals
bool checkTriangleTwist
(
    const bool report,
    const scalar minTwist,
    const polyMesh&,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
);

//- Check for face areas v.s. sum of face-triangle (from face-centre
//  decomposition) areas
bool checkFaceFlatness
(
    const bool report,
    const scalar minFlatness,
    const polyMesh&,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
);

//- Check for small faces
bool checkFaceArea
(
    const bool report,
    const scalar minArea,
    const polyMesh&,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    labelHashSet* setPtr
);

//- Check the area of internal faces v.s. boundary faces
bool checkCellDeterminant
(
    const bool report,
    const scalar minDet,
    const polyMesh&,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    labelHashSet* setPtr
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace polyMeshCheck
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
