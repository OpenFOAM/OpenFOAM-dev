/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "shortEdgeFilter2D.H"

namespace Foam
{
    defineTypeNameAndDebug(shortEdgeFilter2D, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::shortEdgeFilter2D::addRegion
(
    const label regionI,
    DynamicList<label>& bPointRegions
) const
{
    if (bPointRegions.empty())
    {
        bPointRegions.append(regionI);
    }
    else if (findIndex(bPointRegions, regionI) == -1)
    {
        bPointRegions.append(regionI);
    }
}


void Foam::shortEdgeFilter2D::assignBoundaryPointRegions
(
    List<DynamicList<label>>& boundaryPointRegions
) const
{
    forAllConstIter(EdgeMap<label>, mapEdgesRegion_, iter)
    {
        const edge& e = iter.key();
        const label& regionI = iter();

        const label startI = e.start();
        const label endI = e.end();

        addRegion(regionI, boundaryPointRegions[startI]);
        addRegion(regionI, boundaryPointRegions[endI]);
    }
}


void Foam::shortEdgeFilter2D::updateEdgeRegionMap
(
    const MeshedSurface<face>& surfMesh,
    const List<DynamicList<label>>& boundaryPtRegions,
    const labelList& surfPtToBoundaryPt,
    EdgeMap<label>& mapEdgesRegion,
    labelList& patchSizes
) const
{
    EdgeMap<label> newMapEdgesRegion(mapEdgesRegion.size());

    const edgeList& edges = surfMesh.edges();
    const labelList& meshPoints = surfMesh.meshPoints();

    patchSizes.setSize(patchNames_.size(), 0);
    patchSizes = 0;

    forAll(edges, edgeI)
    {
        if (surfMesh.isInternalEdge(edgeI))
        {
            continue;
        }

        const edge& e = edges[edgeI];

        const label startI = meshPoints[e[0]];
        const label endI = meshPoints[e[1]];

        label region = -1;

        const DynamicList<label> startPtRegions =
            boundaryPtRegions[surfPtToBoundaryPt[startI]];
        const DynamicList<label> endPtRegions =
            boundaryPtRegions[surfPtToBoundaryPt[endI]];

        if (startPtRegions.size() > 1 && endPtRegions.size() > 1)
        {
            region = startPtRegions[0];

            WarningInFunction
                << "Both points in edge are in different regions."
                << " Assigning edge to region " << region
                << endl;
        }
        else if (startPtRegions.size() > 1 || endPtRegions.size() > 1)
        {
            region =
            (
                startPtRegions.size() > 1
              ? endPtRegions[0]
              : startPtRegions[0]
            );
        }
        else if
        (
            startPtRegions[0] == endPtRegions[0]
         && startPtRegions[0] != -1
        )
        {
            region = startPtRegions[0];
        }

        if (region != -1)
        {
            newMapEdgesRegion.insert(e, region);
            patchSizes[region]++;
        }
    }

    mapEdgesRegion.transfer(newMapEdgesRegion);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shortEdgeFilter2D::shortEdgeFilter2D
(
    const Foam::CV2D& cv2Dmesh,
    const dictionary& dict
)
:
    shortEdgeFilterFactor_(dict.lookup<scalar>("shortEdgeFilterFactor")),
    edgeAttachedToBoundaryFactor_
    (
        dict.lookupOrDefault<scalar>("edgeAttachedToBoundaryFactor", 2.0)
    ),
    patchNames_(wordList()),
    patchSizes_(labelList()),
    mapEdgesRegion_(),
    indirectPatchEdge_()
{
    point2DField points2D;
    faceList faces;

    cv2Dmesh.calcDual
    (
        points2D,
        faces,
        patchNames_,
        patchSizes_,
        mapEdgesRegion_,
        indirectPatchEdge_
    );

    pointField points(points2D.size());
    forAll(points, ip)
    {
        points[ip] = cv2Dmesh.toPoint3D(points2D[ip]);
    }

    if (debug)
    {
        OFstream str("indirectPatchEdges.obj");
        label count = 0;

        Info<< "Writing indirectPatchEdges to " << str.name() << endl;

        forAllConstIter(EdgeMap<label>, indirectPatchEdge_, iter)
        {
            const edge& e = iter.key();

            meshTools::writeOBJ
            (
                str,
                points[e.start()],
                points[e.end()],
                count
            );
        }
    }

    points2D.clear();

    ms_ = MeshedSurface<face>(move(points), move(faces));

    Info<< "Meshed surface stats before edge filtering :" << endl;
    ms_.writeStats(Info);

    if (debug)
    {
        writeInfo(Info);

        ms_.write("MeshedSurface_preFilter.obj");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shortEdgeFilter2D::~shortEdgeFilter2D()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::shortEdgeFilter2D::filter()
{
    // These are global indices.
    const pointField& points = ms_.points();
    const edgeList& edges = ms_.edges();
    const faceList& faces = ms_.faces();
    const labelList& meshPoints = ms_.meshPoints();
    const labelList& boundaryPoints = ms_.boundaryPoints();

    label maxChain = 0;
    label nPointsToRemove = 0;

    labelList pointsToRemove(ms_.points().size(), -1);

    // List of number of vertices in a face.
    labelList newFaceVertexCount(faces.size(), -1);
    forAll(faces, facei)
    {
        newFaceVertexCount[facei] = faces[facei].size();
    }

    // Check if the point is a boundary point. Flag if it is so that
    // it will not be deleted.
    List<DynamicList<label>> boundaryPointRegions
    (
        points.size(),
        DynamicList<label>()
    );
    assignBoundaryPointRegions(boundaryPointRegions);

    // Check if an edge has a boundary point. It it does the edge length
    // will be doubled when working out its length.
    Info<< "    Marking edges attached to boundaries." << endl;
    boolList edgeAttachedToBoundary(edges.size(), false);
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const label startVertex = e.start();
        const label endVertex = e.end();

        forAll(boundaryPoints, bPoint)
        {
            if
            (
                boundaryPoints[bPoint] == startVertex
             || boundaryPoints[bPoint] == endVertex
            )
            {
                edgeAttachedToBoundary[edgeI] = true;
            }
        }
    }

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        // get the vertices of that edge.
        const label startVertex = e.start();
        const label endVertex = e.end();

        scalar edgeLength =
            mag
            (
                points[meshPoints[startVertex]]
              - points[meshPoints[endVertex]]
            );

        if (edgeAttachedToBoundary[edgeI])
        {
            edgeLength *= edgeAttachedToBoundaryFactor_;
        }

        scalar shortEdgeFilterValue = 0.0;

        const labelList& psEdges = ms_.pointEdges()[startVertex];
        const labelList& peEdges = ms_.pointEdges()[endVertex];

        forAll(psEdges, psEdgeI)
        {
            const edge& psE = edges[psEdges[psEdgeI]];
            if (edgeI != psEdges[psEdgeI])
            {
                shortEdgeFilterValue +=
                    mag
                    (
                        points[meshPoints[psE.start()]]
                      - points[meshPoints[psE.end()]]
                    );
            }
        }

        forAll(peEdges, peEdgeI)
        {
            const edge& peE = edges[peEdges[peEdgeI]];
            if (edgeI != peEdges[peEdgeI])
            {
                shortEdgeFilterValue +=
                    mag
                    (
                        points[meshPoints[peE.start()]]
                      - points[meshPoints[peE.end()]]
                    );
            }
        }

        shortEdgeFilterValue *=
            shortEdgeFilterFactor_
           /(psEdges.size() + peEdges.size() - 2);

        edge lookupInPatchEdge
        (
            meshPoints[startVertex],
            meshPoints[endVertex]
        );

        if
        (
            edgeLength < shortEdgeFilterValue
         || indirectPatchEdge_.found(lookupInPatchEdge)
        )
        {
            bool flagDegenerateFace = false;
            const labelList& pFaces = ms_.pointFaces()[startVertex];

            forAll(pFaces, pFacei)
            {
                const face& f = ms_.localFaces()[pFaces[pFacei]];
                forAll(f, fp)
                {
                    // If the edge is part of this face...
                    if (f[fp] == endVertex)
                    {
                        // If deleting vertex would create a triangle, don't!
                        if (newFaceVertexCount[pFaces[pFacei]] < 4)
                        {
                            flagDegenerateFace = true;
                        }
                        else
                        {
                            newFaceVertexCount[pFaces[pFacei]]--;
                        }
                    }
                    // If the edge is not part of this face...
                    else
                    {
                        // Deleting vertex of a triangle...
                        if (newFaceVertexCount[pFaces[pFacei]] < 3)
                        {
                            flagDegenerateFace = true;
                        }
                    }
                }
            }

            // This if statement determines whether a point should be deleted.
            if
            (
                pointsToRemove[meshPoints[startVertex]] == -1
             && pointsToRemove[meshPoints[endVertex]] == -1
             && !flagDegenerateFace
            )
            {
                const DynamicList<label>& startVertexRegions =
                    boundaryPointRegions[meshPoints[startVertex]];
                const DynamicList<label>& endVertexRegions =
                    boundaryPointRegions[meshPoints[endVertex]];

                if (startVertexRegions.size() && endVertexRegions.size())
                {
                    if (startVertexRegions.size() > 1)
                    {
                        pointsToRemove[meshPoints[endVertex]] =
                            meshPoints[startVertex];
                    }
                    else
                    {
                        pointsToRemove[meshPoints[startVertex]] =
                            meshPoints[endVertex];
                    }
                }
                else if (startVertexRegions.size())
                {
                    pointsToRemove[meshPoints[endVertex]] =
                        meshPoints[startVertex];
                }
                else
                {
                    pointsToRemove[meshPoints[startVertex]] =
                        meshPoints[endVertex];
                }

                ++nPointsToRemove;
            }
        }
    }

    label totalNewPoints = points.size() - nPointsToRemove;

    pointField newPoints(totalNewPoints, Zero);
    labelList newPointNumbers(points.size(), -1);
    label numberRemoved = 0;

    // Maintain addressing from new to old point field
    labelList newPtToOldPt(totalNewPoints, -1);

    forAll(points, pointi)
    {
        // If the point is NOT going to be removed.
        if (pointsToRemove[pointi] == -1)
        {
            newPoints[pointi - numberRemoved] = points[pointi];
            newPointNumbers[pointi] =  pointi - numberRemoved;
            newPtToOldPt[pointi - numberRemoved] = pointi;
        }
        else
        {
            numberRemoved++;
        }
    }

    // Need a new faceList
    faceList newFaces(faces.size());
    label newFacei = 0;

    labelList newFace;
    label newFaceSize = 0;

    // Now need to iterate over the faces and remove points. Global index.
    forAll(faces, facei)
    {
        const face& f = faces[facei];

        newFace.clear();
        newFace.setSize(f.size());
        newFaceSize = 0;

        forAll(f, fp)
        {
            label pointi = f[fp];
            // If not removing the point, then add it to the new face.
            if (pointsToRemove[pointi] == -1)
            {
                newFace[newFaceSize++] = newPointNumbers[pointi];
            }
            else
            {
                label newPointi = pointsToRemove[pointi];
                // Replace deleted point with point that it is being
                // collapsed to.
                if
                (
                    f.nextLabel(fp) != newPointi
                 && f.prevLabel(fp) != newPointi
                )
                {
                    label pChain = newPointi;
                    label totalChain = 0;
                    for (label nChain = 0; nChain <= totalChain; ++nChain)
                    {
                        if (newPointNumbers[pChain] != -1)
                        {
                            newFace[newFaceSize++] = newPointNumbers[pChain];
                            newPointNumbers[pointi] = newPointNumbers[pChain];
                            maxChain = max(totalChain, maxChain);
                        }
                        else
                        {
                            WarningInFunction
                                << "Point " << pChain
                                << " marked for deletion as well as point "
                                << pointi << nl
                                << "    Incrementing maxChain by 1 from "
                                << totalChain << " to " << totalChain + 1
                                << endl;
                            totalChain++;
                        }
                        pChain = pointsToRemove[pChain];
                    }
                }
                else
                {
                    if (newPointNumbers[newPointi] != -1)
                    {
                        newPointNumbers[pointi] = newPointNumbers[newPointi];
                    }
                }
            }
        }

        newFace.setSize(newFaceSize);

        if (newFace.size() > 2)
        {
            newFaces[newFacei++] = face(newFace);
        }
        else
        {
            FatalErrorInFunction
                << "Only " << newFace.size() << " in face " << facei
                << exit(FatalError);
        }
    }

    newFaces.setSize(newFacei);

    MeshedSurface<face> fMesh
    (
        move(newPoints),
        move(newFaces),
        List<surfZone>()
    );

    updateEdgeRegionMap
    (
        fMesh,
        boundaryPointRegions,
        newPtToOldPt,
        mapEdgesRegion_,
        patchSizes_
    );

    forAll(newPointNumbers, pointi)
    {
        if (newPointNumbers[pointi] == -1)
        {
            WarningInFunction
                << pointi << " will be deleted and " << newPointNumbers[pointi]
                << ", so it will not be replaced. "
                << "This will cause edges to be deleted." << endl;
        }
    }

    ms_.transfer(fMesh);

    if (debug)
    {
        Info<< "    Maximum number of chained collapses = " << maxChain << endl;

        writeInfo(Info);
    }
}


void Foam::shortEdgeFilter2D::writeInfo(Ostream& os)
{
    os  << "Short Edge Filtering Information:" << nl
        << "           shortEdgeFilterFactor: " << shortEdgeFilterFactor_ << nl
        << "    edgeAttachedToBoundaryFactor: " << edgeAttachedToBoundaryFactor_
        << endl;

    forAll(patchNames_, patchi)
    {
        os  << "    Patch " << patchNames_[patchi]
            << ", size " << patchSizes_[patchi] << endl;
    }

    os  << "    There are " << mapEdgesRegion_.size()
        << " boundary edges." << endl;

    os  << "    Mesh Info:" << nl
        << "        Points:       " << ms_.nPoints() << nl
        << "        Faces:        " << ms_.size() << nl
        << "        Edges:        " << ms_.nEdges() << nl
        << "            Internal: " << ms_.nInternalEdges() << nl
        << "            External: " << ms_.nEdges() - ms_.nInternalEdges()
        << endl;
}


// ************************************************************************* //
