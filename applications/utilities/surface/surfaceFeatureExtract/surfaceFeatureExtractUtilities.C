/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Extracts and writes surface features to file. All but the basic feature
    extraction is WIP.

\*---------------------------------------------------------------------------*/

#include "surfaceFeatureExtract.H"
#include "Time.H"
#include "tensor2D.H"
#include "symmTensor2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::deleteBox
(
    const triSurface& surf,
    const boundBox& bb,
    const bool removeInside,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    forAll(edgeStat, edgei)
    {
        const point eMid = surf.edges()[edgei].centre(surf.localPoints());

        if (removeInside ? bb.contains(eMid) : !bb.contains(eMid))
        {
            edgeStat[edgei] = surfaceFeatures::NONE;
        }
    }
}


void Foam::deleteEdges
(
    const triSurface& surf,
    const plane& cutPlane,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    const pointField& points = surf.points();
    const labelList& meshPoints = surf.meshPoints();

    forAll(edgeStat, edgei)
    {
        const edge& e = surf.edges()[edgei];
        const point& p0 = points[meshPoints[e.start()]];
        const point& p1 = points[meshPoints[e.end()]];
        const linePointRef line(p0, p1);

        // If edge does not intersect the plane, delete.
        scalar intersect = cutPlane.lineIntersect(line);

        point featPoint = intersect * (p1 - p0) + p0;

        if (!line.insideBoundBox(featPoint))
        {
            edgeStat[edgei] = surfaceFeatures::NONE;
        }
    }
}




void Foam::deleteNonManifoldEdges
(
    const triSurface& surf,
    const scalar tol,
    const scalar includedAngle,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    forAll(edgeStat, edgei)
    {
        const labelList& eFaces = surf.edgeFaces()[edgei];

        if
        (
            eFaces.size() > 2
            && edgeStat[edgei] == surfaceFeatures::REGION
            && (eFaces.size() % 2) == 0
        )
        {
            edgeStat[edgei] = checkNonManifoldEdge
            (
                surf,
                tol,
                includedAngle,
                edgei
            );
        }
    }
}


void Foam::writeStats(const extendedFeatureEdgeMesh& fem, Ostream& os)
{
    os  << "    points : " << fem.points().size() << nl
        << "    of which" << nl
        << "        convex             : "
        << fem.concaveStart() << nl
        << "        concave            : "
        << (fem.mixedStart() - fem.concaveStart()) << nl
        << "        mixed              : "
        << (fem.nonFeatureStart() - fem.mixedStart()) << nl
        << "        non-feature        : "
        << (fem.points().size() - fem.nonFeatureStart()) << nl
        << "    edges  : " << fem.edges().size() << nl
        << "    of which" << nl
        << "        external edges     : "
        << fem.internalStart() << nl
        << "        internal edges     : "
        << (fem.flatStart() - fem.internalStart()) << nl
        << "        flat edges         : "
        << (fem.openStart() - fem.flatStart()) << nl
        << "        open edges         : "
        << (fem.multipleStart() - fem.openStart()) << nl
        << "        multiply connected : "
        << (fem.edges().size() - fem.multipleStart()) << endl;
}


// ************************************************************************* //
