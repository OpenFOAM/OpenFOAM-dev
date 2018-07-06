/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

\*---------------------------------------------------------------------------*/

#include "PatchTools.H"
#include "SortableList.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>

Foam::labelListList
Foam::PatchTools::sortedEdgeFaces
(
    const PrimitivePatch<Face, FaceList, PointField, PointType>& p
)
{
    const edgeList& edges = p.edges();
    const labelListList& edgeFaces = p.edgeFaces();
    const List<Face>& localFaces = p.localFaces();
    const Field<PointType>& localPoints = p.localPoints();

    // create the lists for the various results. (resized on completion)
    labelListList sortedEdgeFaces(edgeFaces.size());

    forAll(edgeFaces, edgeI)
    {
        const labelList& faceNbs = edgeFaces[edgeI];

        if (faceNbs.size() > 2)
        {
            // Get point on edge and normalized direction of edge (= e2 base
            // of our coordinate system)
            const edge& e = edges[edgeI];

            const point& edgePt = localPoints[e.start()];

            vector e2 = e.vec(localPoints);
            e2 /= mag(e2) + vSmall;

            // Get the vertex on 0th face that forms a vector with the first
            // edge point that has the largest angle with the edge
            const Face& f0 = localFaces[faceNbs[0]];

            scalar maxAngle = great;
            vector maxAngleEdgeDir(vector::max);

            forAll(f0, fpI)
            {
                if (f0[fpI] != e.start())
                {
                    vector faceEdgeDir = localPoints[f0[fpI]] - edgePt;
                    faceEdgeDir /= mag(faceEdgeDir) + vSmall;

                    const scalar angle = e2 & faceEdgeDir;

                    if (mag(angle) < maxAngle)
                    {
                        maxAngle = angle;
                        maxAngleEdgeDir = faceEdgeDir;
                    }
                }
            }

            // Get vector normal both to e2 and to edge from opposite vertex
            // to edge (will be x-axis of our coordinate system)
            vector e0 = e2 ^ maxAngleEdgeDir;
            e0 /= mag(e0) + vSmall;

            // Get y-axis of coordinate system
            vector e1 = e2 ^ e0;

            SortableList<scalar> faceAngles(faceNbs.size());

            // e0 is reference so angle is 0
            faceAngles[0] = 0;

            for (label nbI = 1; nbI < faceNbs.size(); nbI++)
            {
                // Get the vertex on face that forms a vector with the first
                // edge point that has the largest angle with the edge
                const Face& f = localFaces[faceNbs[nbI]];

                maxAngle = great;
                maxAngleEdgeDir = vector::max;

                forAll(f, fpI)
                {
                    if (f[fpI] != e.start())
                    {
                        vector faceEdgeDir = localPoints[f[fpI]] - edgePt;
                        faceEdgeDir /= mag(faceEdgeDir) + vSmall;

                        const scalar angle = e2 & faceEdgeDir;

                        if (mag(angle) < maxAngle)
                        {
                            maxAngle = angle;
                            maxAngleEdgeDir = faceEdgeDir;
                        }
                    }
                }

                vector vec = e2 ^ maxAngleEdgeDir;
                vec /= mag(vec) + vSmall;

                faceAngles[nbI] = pseudoAngle
                (
                    e0,
                    e1,
                    vec
                );
            }

            faceAngles.sort();

            sortedEdgeFaces[edgeI] = UIndirectList<label>
            (
                faceNbs,
                faceAngles.indices()
            );
        }
        else
        {
            // No need to sort. Just copy.
            sortedEdgeFaces[edgeI] = faceNbs;
        }
    }

    return sortedEdgeFaces;
}


// ************************************************************************* //
