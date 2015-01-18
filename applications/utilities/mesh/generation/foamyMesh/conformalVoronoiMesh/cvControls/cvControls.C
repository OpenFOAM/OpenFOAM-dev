/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "cvControls.H"
#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cvControls::cvControls
(
    const dictionary& foamyHexMeshDict
)
:
    foamyHexMeshDict_(foamyHexMeshDict)
{
    // Surface conformation controls

    const dictionary& surfDict
    (
        foamyHexMeshDict_.subDict("surfaceConformation")
    );

    pointPairDistanceCoeff_ = readScalar
    (
        surfDict.lookup("pointPairDistanceCoeff")
    );

    mixedFeaturePointPPDistanceCoeff_ = readScalar
    (
        surfDict.lookup("mixedFeaturePointPPDistanceCoeff")
    );

    featurePointExclusionDistanceCoeff_ = readScalar
    (
        surfDict.lookup("featurePointExclusionDistanceCoeff")
    );

    featureEdgeExclusionDistanceCoeff_ = readScalar
    (
        surfDict.lookup("featureEdgeExclusionDistanceCoeff")
    );


    surfaceSearchDistanceCoeff_ = readScalar
    (
        surfDict.lookup("surfaceSearchDistanceCoeff")
    );

    maxSurfaceProtrusionCoeff_ = readScalar
    (
        surfDict.lookup("maxSurfaceProtrusionCoeff")
    );

    maxQuadAngle_ = readScalar(surfDict.lookup("maxQuadAngle"));

    surfaceConformationRebuildFrequency_ = max
    (
        1,
        readLabel(surfDict.lookup("surfaceConformationRebuildFrequency"))
    );


    const dictionary& featurePointControlsDict
    (
        surfDict.subDict("featurePointControls")
    );

    specialiseFeaturePoints_ = Switch
    (
        featurePointControlsDict.lookup("specialiseFeaturePoints")
    );

    guardFeaturePoints_ = Switch
    (
        featurePointControlsDict.lookup("guardFeaturePoints")
    );

    edgeAiming_ = Switch
    (
        featurePointControlsDict.lookup("edgeAiming")
    );

    if (!guardFeaturePoints_)
    {
        snapFeaturePoints_ = Switch
        (
            featurePointControlsDict.lookup("snapFeaturePoints")
        );
    }

    circulateEdges_ = Switch
    (
        featurePointControlsDict.lookup("circulateEdges")
    );

    // Controls for coarse surface conformation

    const dictionary& conformationControlsDict
    (
        surfDict.subDict("conformationControls")
    );

    surfacePtExclusionDistanceCoeff_ = readScalar
    (
        conformationControlsDict.lookup("surfacePtExclusionDistanceCoeff")
    );

    edgeSearchDistCoeffSqr_ = sqr
    (
        readScalar
        (
            conformationControlsDict.lookup("edgeSearchDistCoeff")
        )
    );

    surfacePtReplaceDistCoeffSqr_ = sqr
    (
        readScalar
        (
            conformationControlsDict.lookup("surfacePtReplaceDistCoeff")
        )
    );

    maxConformationIterations_ = readLabel
    (
        conformationControlsDict.lookup("maxIterations")
    );

    iterationToInitialHitRatioLimit_ = readScalar
    (
        conformationControlsDict.lookup("iterationToInitialHitRatioLimit")
    );


    // Motion control controls

    const dictionary& motionDict(foamyHexMeshDict_.subDict("motionControl"));

    defaultCellSize_ = readScalar(motionDict.lookup("defaultCellSize"));

    minimumCellSize_ =
        readScalar(motionDict.lookup("minimumCellSizeCoeff"))*defaultCellSize_;

    objOutput_ = Switch(motionDict.lookupOrDefault<Switch>("objOutput", false));

    timeChecks_ = Switch
    (
        motionDict.lookupOrDefault<Switch>("timeChecks", false)
    );

    printVertexInfo_ = Switch
    (
        motionDict.lookupOrDefault<Switch>("printVertexInfo", false)
    );

    if (Pstream::parRun())
    {
        maxLoadUnbalance_ = readScalar(motionDict.lookup("maxLoadUnbalance"));
    }
    else
    {
        maxLoadUnbalance_ = -1;
    }

    cosAlignmentAcceptanceAngle_ = cos
    (
        degToRad(readScalar(motionDict.lookup("alignmentAcceptanceAngle")))
    );


    // Point removal criteria

    const dictionary& insertionDict
    (
        motionDict.subDict("pointInsertionCriteria")
    );

    insertionDistCoeff_ = readScalar
    (
        insertionDict.lookup("cellCentreDistCoeff")
    );

    faceAreaRatioCoeff_ = readScalar
    (
        insertionDict.lookup("faceAreaRatioCoeff")
    );

    cosInsertionAcceptanceAngle_ = cos
    (
        degToRad(readScalar(insertionDict.lookup("acceptanceAngle")))
    );

    // Point removal criteria

    const dictionary& removalDict
    (
        motionDict.subDict("pointRemovalCriteria")
    );

    removalDistCoeff_ = readScalar
    (
        removalDict.lookup("cellCentreDistCoeff")
    );

    // polyMesh filtering controls

    const dictionary& filteringDict
    (
        foamyHexMeshDict_.subDict("polyMeshFiltering")
    );

    filterEdges_ = Switch
    (
        filteringDict.lookupOrDefault<Switch>("filterEdges", true)
    );

    filterFaces_ = Switch
    (
        filteringDict.lookupOrDefault<Switch>("filterFaces", false)
    );

    if (filterFaces_)
    {
        filterEdges_ = Switch::ON;
    }

    writeTetDualMesh_ = Switch(filteringDict.lookup("writeTetDualMesh"));

    writeCellShapeControlMesh_ =
        Switch(filteringDict.lookup("writeCellShapeControlMesh"));

    if (Pstream::parRun())
    {
        writeBackgroundMeshDecomposition_ =
            Switch(filteringDict.lookup("writeBackgroundMeshDecomposition"));
    }
    else
    {
        writeBackgroundMeshDecomposition_ = Switch(false);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cvControls::~cvControls()
{}


// ************************************************************************* //
