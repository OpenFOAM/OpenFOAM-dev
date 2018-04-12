/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

Application
    surfaceFeatureExtract

Description
    Utility functions for surfaceFeatureExtract

\*---------------------------------------------------------------------------*/

#include "surfaceFeatures.H"
#include "extendedFeatureEdgeMesh.H"
#include "triSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "pointIndexHitList.H"
#include "plane.H"
#include "triadField.H"

namespace Foam
{
    extern const scalar internalAngleTolerance;
    extern const scalar internalToleranceCosAngle;

    extern const scalar externalAngleTolerance;
    extern const scalar externalToleranceCosAngle;

    //- Deletes all edges inside/outside bounding box from set.
    void deleteBox
    (
        const triSurface& surf,
        const boundBox& bb,
        const bool removeInside,
        List<surfaceFeatures::edgeStatus>& edgeStat
    );

    //- Deletes all edges inside/outside bounding box from set.
    void deleteEdges
    (
        const triSurface& surf,
        const plane& cutPlane,
        List<surfaceFeatures::edgeStatus>& edgeStat
    );

    //- Divide into multiple normal bins
    //  - return REGION if != 2 normals
    //  - return REGION if 2 normals that make feature angle
    //  - otherwise return NONE and set normals,bins
    surfaceFeatures::edgeStatus checkNonManifoldEdge
    (
        const triSurface& surf,
        const scalar tol,
        const scalar includedAngle,
        const label edgei
    );

    void deleteNonManifoldEdges
    (
        const triSurface& surf,
        const scalar tol,
        const scalar includedAngle,
        List<surfaceFeatures::edgeStatus>& edgeStat
    );

    void writeStats(const extendedFeatureEdgeMesh& fem, Ostream& os);


    void drawHitProblem
    (
        const label fi,
        const triSurface& surf,
        const point& start,
        const point& p,
        const point& end,
        const pointIndexHitList& hitInfo
    );

    void processHit
    (
        scalar& internalCloseness,
        scalar& externalCloseness,
        const label fi,
        const triSurface& surf,
        const point& start,
        const point& p,
        const point& end,
        const vector& normal,
        const vectorField& normals,
        const pointIndexHitList& hitInfo
    );

    Pair<tmp<triSurfaceScalarField>> extractCloseness
    (
        const triSurfaceMesh& surf
    );

    Pair<tmp<triSurfacePointScalarField>> extractPointCloseness
    (
        const triSurfaceMesh& surf
    );
}


// ************************************************************************* //
