/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::dictionary& Foam::cvControls::foamyHexMeshDict() const
{
    return foamyHexMeshDict_;
}


inline Foam::scalar Foam::cvControls::pointPairDistanceCoeff() const
{
    return pointPairDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::mixedFeaturePointPPDistanceCoeff() const
{
    return mixedFeaturePointPPDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::featurePointExclusionDistanceCoeff() const
{
    return featurePointExclusionDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::featureEdgeExclusionDistanceCoeff() const
{
    return featureEdgeExclusionDistanceCoeff_;
}


inline Foam::Switch Foam::cvControls::guardFeaturePoints() const
{
    return guardFeaturePoints_;
}

inline Foam::Switch Foam::cvControls::edgeAiming() const
{
    return edgeAiming_;
}

inline Foam::Switch Foam::cvControls::snapFeaturePoints() const
{
    return snapFeaturePoints_;
}

inline Foam::Switch Foam::cvControls::circulateEdges() const
{
    return circulateEdges_;
}


inline Foam::scalar Foam::cvControls::surfacePtExclusionDistanceCoeff() const
{
    return surfacePtExclusionDistanceCoeff_;
}

inline Foam::Switch Foam::cvControls::specialiseFeaturePoints() const
{
    return specialiseFeaturePoints_;
}

inline Foam::scalar Foam::cvControls::surfaceSearchDistanceCoeff() const
{
    return surfaceSearchDistanceCoeff_;
}


inline Foam::scalar Foam::cvControls::maxSurfaceProtrusionCoeff() const
{
    return maxSurfaceProtrusionCoeff_;
}


inline Foam::scalar Foam::cvControls::maxQuadAngle() const
{
    return maxQuadAngle_;
}


inline Foam::label Foam::cvControls::surfaceConformationRebuildFrequency() const
{
    return surfaceConformationRebuildFrequency_;
}


inline Foam::scalar Foam::cvControls::edgeSearchDistCoeffSqr() const
{
    return edgeSearchDistCoeffSqr_;
}


inline Foam::scalar Foam::cvControls::surfacePtReplaceDistCoeffSqr() const
{
    return surfacePtReplaceDistCoeffSqr_;
}


inline Foam::label Foam::cvControls::maxConformationIterations() const
{
    return maxConformationIterations_;
}


inline Foam::scalar Foam::cvControls::iterationToInitialHitRatioLimit() const
{
    return iterationToInitialHitRatioLimit_;
}


inline Foam::scalar Foam::cvControls::defaultCellSize() const
{
    return defaultCellSize_;
}


inline Foam::scalar Foam::cvControls::minimumCellSize() const
{
    return minimumCellSize_;
}


inline Foam::Switch Foam::cvControls::objOutput() const
{
    return objOutput_;
}


inline Foam::Switch Foam::cvControls::timeChecks() const
{
    return timeChecks_;
}


inline Foam::Switch Foam::cvControls::printVertexInfo() const
{
    return printVertexInfo_;
}


inline Foam::scalar Foam::cvControls::maxLoadUnbalance() const
{
    return maxLoadUnbalance_;
}


inline Foam::scalar Foam::cvControls::cosAlignmentAcceptanceAngle() const
{
    return cosAlignmentAcceptanceAngle_;
}


inline Foam::scalar Foam::cvControls::insertionDistCoeff() const
{
    return insertionDistCoeff_;
}


inline Foam::scalar Foam::cvControls::faceAreaRatioCoeff() const
{
    return faceAreaRatioCoeff_;
}


inline Foam::scalar Foam::cvControls::cosInsertionAcceptanceAngle() const
{
    return cosInsertionAcceptanceAngle_;
}


inline Foam::scalar Foam::cvControls::removalDistCoeff() const
{
    return removalDistCoeff_;
}


inline Foam::Switch Foam::cvControls::filterEdges() const
{
    return filterEdges_;
}

inline Foam::Switch Foam::cvControls::filterFaces() const
{
    return filterFaces_;
}

inline Foam::Switch Foam::cvControls::writeTetDualMesh() const
{
    return writeTetDualMesh_;
}

inline Foam::Switch Foam::cvControls::writeCellShapeControlMesh() const
{
    return writeCellShapeControlMesh_;
}

inline Foam::Switch Foam::cvControls::writeBackgroundMeshDecomposition() const
{
    return writeBackgroundMeshDecomposition_;
}


// ************************************************************************* //
