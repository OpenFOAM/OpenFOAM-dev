/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

    pointPairDistanceCoeff_ =
        surfDict.lookup<scalar>("pointPairDistanceCoeff");

    mixedFeaturePointPPDistanceCoeff_ =
        surfDict.lookup<scalar>("mixedFeaturePointPPDistanceCoeff");

    featurePointExclusionDistanceCoeff_ =
        surfDict.lookup<scalar>("featurePointExclusionDistanceCoeff");

    featureEdgeExclusionDistanceCoeff_ =
        surfDict.lookup<scalar>("featureEdgeExclusionDistanceCoeff");

    surfaceSearchDistanceCoeff_ =
        surfDict.lookup<scalar>("surfaceSearchDistanceCoeff");

    maxSurfaceProtrusionCoeff_ =
        surfDict.lookup<scalar>("maxSurfaceProtrusionCoeff");

    maxQuadAngle_ = surfDict.lookup<scalar>("maxQuadAngle");

    surfaceConformationRebuildFrequency_ = max
    (
        1,
        surfDict.lookup<label>("surfaceConformationRebuildFrequency")
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

    surfacePtExclusionDistanceCoeff_ =
        conformationControlsDict.lookup<scalar>
        (
            "surfacePtExclusionDistanceCoeff"
        );

    edgeSearchDistCoeffSqr_ = sqr
    (
        conformationControlsDict.lookup<scalar>("edgeSearchDistCoeff")
    );

    surfacePtReplaceDistCoeffSqr_ = sqr
    (
        conformationControlsDict.lookup<scalar>("surfacePtReplaceDistCoeff")
    );

    maxConformationIterations_ =
        conformationControlsDict.lookup<label>("maxIterations");

    iterationToInitialHitRatioLimit_ =
        conformationControlsDict.lookup<scalar>
        (
            "iterationToInitialHitRatioLimit"
        );


    // Motion control controls

    const dictionary& motionDict(foamyHexMeshDict_.subDict("motionControl"));

    defaultCellSize_ = motionDict.lookup<scalar>("defaultCellSize");

    minimumCellSize_ =
        motionDict.lookup<scalar>("minimumCellSizeCoeff")*defaultCellSize_;

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
        maxLoadUnbalance_ = motionDict.lookup<scalar>("maxLoadUnbalance");
    }
    else
    {
        maxLoadUnbalance_ = -1;
    }

    cosAlignmentAcceptanceAngle_ = cos
    (
        degToRad(motionDict.lookup<scalar>("alignmentAcceptanceAngle"))
    );


    // Point removal criteria

    const dictionary& insertionDict
    (
        motionDict.subDict("pointInsertionCriteria")
    );

    insertionDistCoeff_ =
        insertionDict.lookup<scalar>("cellCentreDistCoeff");

    faceAreaRatioCoeff_ =
        insertionDict.lookup<scalar>("faceAreaRatioCoeff");

    cosInsertionAcceptanceAngle_ = cos
    (
        degToRad(insertionDict.lookup<scalar>("acceptanceAngle"))
    );

    // Point removal criteria

    const dictionary& removalDict
    (
        motionDict.subDict("pointRemovalCriteria")
    );

    removalDistCoeff_ =
        removalDict.lookup<scalar>("cellCentreDistCoeff");

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
        filterEdges_ = Switch::switchType::on;
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
