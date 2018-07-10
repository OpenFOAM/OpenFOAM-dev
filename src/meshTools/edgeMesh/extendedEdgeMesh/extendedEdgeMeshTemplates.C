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

#include "extendedEdgeMesh.H"
#include "ListListOps.H"
#include "unitConversion.H"
#include "PackedBoolList.H"
#include "PatchTools.H"
#include "searchableBox.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Patch>
void Foam::extendedEdgeMesh::sortPointsAndEdges
(
    const Patch& surf,
    const labelList& featureEdges,
    const labelList& regionFeatureEdges,// subset of featureEdges: inter-region
    const labelList& featurePoints
)
{
    const pointField& sFeatLocalPts(surf.localPoints());
    const edgeList& sFeatEds(surf.edges());
    const labelListList edgeFaces = PatchTools::sortedEdgeFaces(surf);
    const vectorField& faceNormals = surf.faceNormals();
    const labelListList pointEdges = PatchTools::sortedPointEdges(surf);

    // Extract and reorder the data from surfaceFeatures

    // References to the surfaceFeatures data

    // Filling the extendedEdgeMesh with the raw geometrical data.

    label nFeatEds = featureEdges.size();
    label nFeatPts = featurePoints.size();

    DynamicList<point> tmpPts;
    edgeList eds(nFeatEds);
    DynamicList<vector> norms;
    vectorField edgeDirections(nFeatEds);
    labelListList edgeNormals(nFeatEds);
    labelListList normalDirections(nFeatEds);
    DynamicList<label> regionEdges;

    // Keep track of the ordered feature point feature edges
    labelListList featurePointFeatureEdges(nFeatPts);
    forAll(featurePointFeatureEdges, pI)
    {
        featurePointFeatureEdges[pI] = pointEdges[featurePoints[pI]];
    }

    // Mapping between old and new indices, there is entry in the map for each
    // of surf.localPoints, -1 means that this point hasn't been used (yet),
    // >= 0 corresponds to the index
    labelList pointMap(sFeatLocalPts.size(), -1);

    // Mapping between surface edge index and its feature edge index. -1 if it
    // is not a feature edge
    labelList edgeMap(sFeatEds.size(), -1);

    // Noting when the normal of a face has been used so not to duplicate
    labelList faceMap(surf.size(), -1);

    // Collecting the status of edge for subsequent sorting
    List<edgeStatus> edStatus(nFeatEds, NONE);

    forAll(featurePoints, i)
    {
        label sFPI = featurePoints[i];

        tmpPts.append(sFeatLocalPts[sFPI]);

        pointMap[sFPI] = tmpPts.size() - 1;
    }

    // All feature points have been added
    nonFeatureStart_ = tmpPts.size();

    PackedBoolList isRegionFeatureEdge(regionFeatureEdges);

    forAll(featureEdges, i)
    {
        label sFEI = featureEdges[i];

        edgeMap[sFEI] = i;

        const edge& fE = sFeatEds[sFEI];

        edgeDirections[i] = fE.vec(sFeatLocalPts);

        // Check to see if the points have been already used
        if (pointMap[fE.start()] == -1)
        {
             tmpPts.append(sFeatLocalPts[fE.start()]);

             pointMap[fE.start()] = tmpPts.size() - 1;
        }

        eds[i].start() = pointMap[fE.start()];

        if (pointMap[fE.end()] == -1)
        {
            tmpPts.append(sFeatLocalPts[fE.end()]);

            pointMap[fE.end()] = tmpPts.size() - 1;
        }

        eds[i].end() = pointMap[fE.end()];

        // Pick up the faces adjacent to the feature edge
        const labelList& eFaces = edgeFaces[sFEI];

        edgeNormals[i].setSize(eFaces.size());
        normalDirections[i].setSize(eFaces.size());

        forAll(eFaces, j)
        {
            label eFI = eFaces[j];

            // Check to see if the points have been already used
            if (faceMap[eFI] == -1)
            {
                norms.append(faceNormals[eFI]);

                faceMap[eFI] = norms.size() - 1;
            }

            edgeNormals[i][j] = faceMap[eFI];

            const vector cross = (faceNormals[eFI] ^ edgeDirections[i]);
            const vector fC0tofE0 =
                surf[eFI].centre(surf.points())
              - sFeatLocalPts[fE.start()];

            normalDirections[i][j] =
                (
                    (
                        (cross/(mag(cross) + vSmall))
                      & (fC0tofE0/(mag(fC0tofE0)+ vSmall))
                    )
                  > 0.0
                    ? 1
                    : -1
                );
        }

        vector fC0tofC1(Zero);

        if (eFaces.size() == 2)
        {
            fC0tofC1 =
                surf[eFaces[1]].centre(surf.points())
              - surf[eFaces[0]].centre(surf.points());
        }

        edStatus[i] = classifyEdge(norms, edgeNormals[i], fC0tofC1);

        if (isRegionFeatureEdge[i])
        {
            regionEdges.append(i);
        }
    }

    // Populate feature point feature edges
    DynamicList<label> newFeatureEdges;

    forAll(featurePointFeatureEdges, pI)
    {
        const labelList& fpfe = featurePointFeatureEdges[pI];

        newFeatureEdges.setCapacity(fpfe.size());

        forAll(fpfe, eI)
        {
            const label oldEdgeIndex = fpfe[eI];

            const label newFeatureEdgeIndex = edgeMap[oldEdgeIndex];

            if (newFeatureEdgeIndex != -1)
            {
                newFeatureEdges.append(newFeatureEdgeIndex);
            }
        }

        featurePointFeatureEdges[pI].transfer(newFeatureEdges);
    }

    // Reorder the edges by classification
    List<DynamicList<label>> allEds(nEdgeTypes);

    DynamicList<label>& externalEds(allEds[0]);
    DynamicList<label>& internalEds(allEds[1]);
    DynamicList<label>& flatEds(allEds[2]);
    DynamicList<label>& openEds(allEds[3]);
    DynamicList<label>& multipleEds(allEds[4]);

    forAll(eds, i)
    {
        edgeStatus eStat = edStatus[i];

        if (eStat == EXTERNAL)
        {
            externalEds.append(i);
        }
        else if (eStat == INTERNAL)
        {
            internalEds.append(i);
        }
        else if (eStat == FLAT)
        {
            flatEds.append(i);
        }
        else if (eStat == OPEN)
        {
            openEds.append(i);
        }
        else if (eStat == MULTIPLE)
        {
            multipleEds.append(i);
        }
        else if (eStat == NONE)
        {
            FatalErrorInFunction
                << nl << "classifyEdge returned NONE on edge "
                << eds[i]
                << ". There is a problem with definition of this edge."
                << nl << abort(FatalError);
        }
    }

    internalStart_ = externalEds.size();
    flatStart_ = internalStart_ + internalEds.size();
    openStart_ = flatStart_ + flatEds.size();
    multipleStart_ = openStart_ + openEds.size();

    labelList edMap
    (
        ListListOps::combine<labelList>
        (
            allEds,
            accessOp<labelList>()
        )
    );

    edMap = invert(edMap.size(), edMap);

    inplaceReorder(edMap, eds);
    inplaceReorder(edMap, edStatus);
    inplaceReorder(edMap, edgeDirections);
    inplaceReorder(edMap, edgeNormals);
    inplaceReorder(edMap, normalDirections);
    inplaceRenumber(edMap, regionEdges);

    forAll(featurePointFeatureEdges, pI)
    {
        inplaceRenumber(edMap, featurePointFeatureEdges[pI]);
    }

    pointField pts(tmpPts);

    // Initialise the edgeMesh
    edgeMesh::operator=(edgeMesh(pts, eds));

    // Initialise sorted edge related data
    edgeDirections_ = edgeDirections/(mag(edgeDirections) + vSmall);
    edgeNormals_ = edgeNormals;
    normalDirections_ = normalDirections;
    regionEdges_ = regionEdges;

    // Normals are all now found and indirectly addressed, can also be stored
    normals_ = vectorField(norms);


    // Reorder the feature points by classification

    List<DynamicList<label>> allPts(3);

    DynamicList<label>& convexPts(allPts[0]);
    DynamicList<label>& concavePts(allPts[1]);
    DynamicList<label>& mixedPts(allPts[2]);

    for (label i = 0; i < nonFeatureStart_; i++)
    {
        pointStatus ptStatus = classifyFeaturePoint(i);

        if (ptStatus == CONVEX)
        {
            convexPts.append(i);
        }
        else if (ptStatus == CONCAVE)
        {
            concavePts.append(i);
        }
        else if (ptStatus == MIXED)
        {
            mixedPts.append(i);
        }
        else if (ptStatus == NONFEATURE)
        {
            FatalErrorInFunction
                << nl << "classifyFeaturePoint returned NONFEATURE on point at "
                << points()[i]
                << ". There is a problem with definition of this feature point."
                << nl << abort(FatalError);
        }
    }

    concaveStart_ = convexPts.size();
    mixedStart_ = concaveStart_ + concavePts.size();

    labelList ftPtMap
    (
        ListListOps::combine<labelList>
        (
            allPts,
            accessOp<labelList>()
        )
    );

    ftPtMap = invert(ftPtMap.size(), ftPtMap);

    // Creating the ptMap from the ftPtMap with identity values up to the size
    // of pts to create an oldToNew map for inplaceReorder

    labelList ptMap(identity(pts.size()));

    forAll(ftPtMap, i)
    {
        ptMap[i] = ftPtMap[i];
    }

    inplaceReorder(ptMap, pts);
    inplaceReorder(ptMap, featurePointFeatureEdges);

    forAll(eds, i)
    {
        inplaceRenumber(ptMap, eds[i]);
    }

    // Reinitialise the edgeMesh with sorted feature points and
    // renumbered edges
    reset(xferMove(pts), xferMove(eds));

    // Generate the featurePointNormals

    labelListList featurePointNormals(nonFeatureStart_);

    for (label i = 0; i < nonFeatureStart_; i++)
    {
        DynamicList<label> tmpFtPtNorms;

        const labelList& ptEds = edgeMesh::pointEdges()[i];

        forAll(ptEds, j)
        {
            const labelList& ptEdNorms(edgeNormals[ptEds[j]]);

            forAll(ptEdNorms, k)
            {
                if (findIndex(tmpFtPtNorms, ptEdNorms[k]) == -1)
                {
                    bool addNormal = true;

                    // Check that the normal direction is unique at this feature
                    forAll(tmpFtPtNorms, q)
                    {
                        if
                        (
                            (normals_[ptEdNorms[k]] & normals_[tmpFtPtNorms[q]])
                          > cosNormalAngleTol_
                        )
                        {
                            // Parallel to an existing normal, do not add
                            addNormal = false;

                            break;
                        }
                    }

                    if (addNormal)
                    {
                        tmpFtPtNorms.append(ptEdNorms[k]);
                    }
                }
            }
        }

        featurePointNormals[i] = tmpFtPtNorms;
    }

    featurePointNormals_ = featurePointNormals;
    featurePointEdges_ = featurePointFeatureEdges;
}


// ************************************************************************* //
