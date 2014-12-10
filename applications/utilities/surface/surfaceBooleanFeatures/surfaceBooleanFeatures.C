/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    surfaceBooleanFeatures

Description

    Generates the extendedFeatureEdgeMesh for the interface between a boolean
    operation on two surfaces.  Assumes that the orientation of the surfaces is
    correct:

    + if the operation is union or intersection, that both surface's normals
      (n) have the same orientation with respect to a point, i.e. surfaces and b
      are orientated the same with respect to point x:

    @verbatim
       _______
      |       |--> n
      |    ___|___             x
      |a  |   |   |--> n
      |___|___|  b|
          |       |
          |_______|

    @endverbatim

    + if the operation is a subtraction, the surfaces should be oppositely
    oriented with respect to a point, i.e. for (a - b), then b's orientation
    should be such that x is "inside", and a's orientation such that x is
    "outside"

    @verbatim
       _______
      |       |--> n
      |    ___|___             x
      |a  |   |   |
      |___|___|  b|
          |  n <--|
          |_______|

    @endverbatim

    When the operation is peformed - for union, all of the edges generates where
    one surfaces cuts another are all "internal" for union, and "external" for
    intersection, b - a and a - b.  This has been assumed, formal (dis)proof is
    invited.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "Time.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "triSurfaceSearch.H"
#include "OFstream.H"
#include "booleanSurface.H"
#include "edgeIntersections.H"
#include "meshTools.H"
#include "labelPair.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool intersectSurfaces
(
    triSurface& surf,
    edgeIntersections& edgeCuts
)
{
    bool hasMoved = false;

    for (label iter = 0; iter < 10; iter++)
    {
        Info<< "Determining intersections of surface edges with itself" << endl;

        // Determine surface edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters = 0;
        {
            triSurfaceSearch querySurf(surf);

            scalarField surfPointTol
            (
                1e-3*edgeIntersections::minEdgeLength(surf)
            );

            // Determine raw intersections
            edgeCuts = edgeIntersections
            (
                surf,
                querySurf,
                surfPointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points(surf.points());

                nIters =
                    edgeCuts.removeDegenerates
                    (
                        5,              // max iterations
                        surf,
                        querySurf,
                        surfPointTol,
                        points         // work array
                    );

                if (nIters != 0)
                {
                    // Update geometric quantities
                    surf.movePoints(points);
                    hasMoved = true;
                }
            }
        }
    }

    if (hasMoved)
    {
        fileName newFile("surf.obj");
        Info<< "Surface has been moved. Writing to " << newFile << endl;
        surf.write(newFile);
    }

    return hasMoved;
}


// Keep on shuffling surface points until no more degenerate intersections.
// Moves both surfaces and updates set of edge cuts.
bool intersectSurfaces
(
    triSurface& surf1,
    edgeIntersections& edgeCuts1,
    triSurface& surf2,
    edgeIntersections& edgeCuts2
)
{
    bool hasMoved1 = false;
    bool hasMoved2 = false;

    for (label iter = 0; iter < 10; iter++)
    {
        Info<< "Determining intersections of surf1 edges with surf2"
            << " faces" << endl;

        // Determine surface1 edge intersections. Allow surface to be moved.

        // Number of iterations needed to resolve degenerates
        label nIters1 = 0;
        {
            triSurfaceSearch querySurf2(surf2);

            scalarField surf1PointTol
            (
                1e-3*edgeIntersections::minEdgeLength(surf1)
            );

            // Determine raw intersections
            edgeCuts1 = edgeIntersections
            (
                surf1,
                querySurf2,
                surf1PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points1(surf1.points());

                nIters1 =
                    edgeCuts1.removeDegenerates
                    (
                        5,              // max iterations
                        surf1,
                        querySurf2,
                        surf1PointTol,
                        points1         // work array
                    );

                if (nIters1 != 0)
                {
                    // Update geometric quantities
                    surf1.movePoints(points1);
                    hasMoved1 = true;
                }
            }
        }

        Info<< "Determining intersections of surf2 edges with surf1"
            << " faces" << endl;

        label nIters2 = 0;
        {
            triSurfaceSearch querySurf1(surf1);

            scalarField surf2PointTol
            (
                1e-3*edgeIntersections::minEdgeLength(surf2)
            );

            // Determine raw intersections
            edgeCuts2 = edgeIntersections
            (
                surf2,
                querySurf1,
                surf2PointTol
            );

            // Shuffle a bit to resolve degenerate edge-face hits
            {
                pointField points2(surf2.points());

                nIters2 =
                    edgeCuts2.removeDegenerates
                    (
                        5,              // max iterations
                        surf2,
                        querySurf1,
                        surf2PointTol,
                        points2         // work array
                    );

                if (nIters2 != 0)
                {
                    // Update geometric quantities
                    surf2.movePoints(points2);
                    hasMoved2 = true;
                }
            }
        }

        if (nIters1 == 0 && nIters2 == 0)
        {
            Info<< "** Resolved all intersections to be proper edge-face pierce"
                << endl;
            break;
        }
    }

    if (hasMoved1)
    {
        fileName newFile("surf1.obj");
        Info<< "Surface 1 has been moved. Writing to " << newFile
            << endl;
        surf1.write(newFile);
    }

    if (hasMoved2)
    {
        fileName newFile("surf2.obj");
        Info<< "Surface 2 has been moved. Writing to " << newFile
            << endl;
        surf2.write(newFile);
    }

    return hasMoved1 || hasMoved2;
}


label calcNormalDirection
(
    const vector& normal,
    const vector& otherNormal,
    const vector& edgeDir,
    const vector& faceCentre,
    const vector& pointOnEdge
)
{
    vector cross = (normal ^ edgeDir);
    cross /= mag(cross);

    vector fC0tofE0 = faceCentre - pointOnEdge;
    fC0tofE0 /= mag(fC0tofE0);

    label nDir = ((cross & fC0tofE0) > 0.0 ? 1 : -1);

    nDir *= ((otherNormal & fC0tofE0) > 0.0 ? -1 : 1);

    return nDir;
}


void calcEdgeCuts
(
    triSurface& surf1,
    triSurface& surf2,
    const bool perturb,
    edgeIntersections& edge1Cuts,
    edgeIntersections& edge2Cuts
)
{
    if (perturb)
    {
        intersectSurfaces
        (
            surf1,
            edge1Cuts,
            surf2,
            edge2Cuts
        );
    }
    else
    {
        triSurfaceSearch querySurf2(surf2);

        Info<< "Determining intersections of surf1 edges with surf2 faces"
            << endl;

        edge1Cuts = edgeIntersections
        (
            surf1,
            querySurf2,
            1e-3*edgeIntersections::minEdgeLength(surf1)
        );

        triSurfaceSearch querySurf1(surf1);

        Info<< "Determining intersections of surf2 edges with surf1 faces"
            << endl;

        edge2Cuts = edgeIntersections
        (
            surf2,
            querySurf1,
            1e-3*edgeIntersections::minEdgeLength(surf2)
        );
    }
}


void calcFeaturePoints(const pointField& points, const edgeList& edges)
{
    edgeMesh eMesh(points, edges);

    const labelListList& pointEdges = eMesh.pointEdges();


    // Get total number of feature points
    label nFeaturePoints = 0;
    forAll(pointEdges, pI)
    {
        const labelList& pEdges = pointEdges[pI];

        if (pEdges.size() == 1)
        {
            nFeaturePoints++;
        }
    }


    // Calculate addressing from feature point to cut point and cut edge
    labelList featurePointToCutPoint(nFeaturePoints);
    labelList featurePointToCutEdge(nFeaturePoints);

    label nFeatPts = 0;
    forAll(pointEdges, pI)
    {
        const labelList& pEdges = pointEdges[pI];

        if (pEdges.size() == 1)
        {
            featurePointToCutPoint[nFeatPts] = pI;
            featurePointToCutEdge[nFeatPts] = pEdges[0];
            nFeatPts++;
        }
    }



    label concaveStart = 0;
    label mixedStart = 0;
    label nonFeatureStart = nFeaturePoints;


    labelListList featurePointNormals(nFeaturePoints);
    labelListList featurePointEdges(nFeaturePoints);
    labelList regionEdges;


}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("action");
    argList::validArgs.append("surface file");
    argList::validArgs.append("surface file");

    argList::addBoolOption
    (
        "surf1Baffle",
        "Mark surface 1 as a baffle"
    );

    argList::addBoolOption
    (
        "surf2Baffle",
        "Mark surface 2 as a baffle"
    );

    argList::addBoolOption
    (
        "perturb",
        "Perturb surface points to escape degenerate intersections"
    );

    argList::addBoolOption
    (
        "invertedSpace",
        "do the surfaces have inverted space orientation, "
        "i.e. a point at infinity is considered inside. "
        "This is only sensible for union and intersection."
    );

    #   include "setRootCase.H"
    #   include "createTime.H"

    word action(args.args()[1]);

    HashTable<booleanSurface::booleanOpType> validActions;
    validActions.insert("intersection", booleanSurface::INTERSECTION);
    validActions.insert("union", booleanSurface::UNION);
    validActions.insert("difference", booleanSurface::DIFFERENCE);

    if (!validActions.found(action))
    {
        FatalErrorIn(args.executable())
            << "Unsupported action " << action << endl
            << "Supported actions:" << validActions.toc() << abort(FatalError);
    }

    fileName surf1Name(args.args()[2]);
    Info<< "Reading surface " << surf1Name << endl;
    triSurface surf1(surf1Name);

    Info<< surf1Name << " statistics:" << endl;
    surf1.writeStats(Info);
    Info<< endl;

    fileName surf2Name(args[3]);
    Info<< "Reading surface " << surf2Name << endl;
    triSurface surf2(surf2Name);

    Info<< surf2Name << " statistics:" << endl;
    surf2.writeStats(Info);
    Info<< endl;

    const bool surf1Baffle = args.optionFound("surf1Baffle");
    const bool surf2Baffle = args.optionFound("surf2Baffle");

    edgeIntersections edge1Cuts;
    edgeIntersections edge2Cuts;

    bool invertedSpace = args.optionFound("invertedSpace");

    if (invertedSpace && validActions[action] == booleanSurface::DIFFERENCE)
    {
        FatalErrorIn(args.executable())
            << "Inverted space only makes sense for union or intersection."
            << exit(FatalError);
    }

    // Calculate the points where the edges are cut by the other surface
    calcEdgeCuts
    (
        surf1,
        surf2,
        args.optionFound("perturb"),
        edge1Cuts,
        edge2Cuts
    );

    // Determine intersection edges from the edge cuts
    surfaceIntersection inter
    (
        surf1,
        edge1Cuts,
        surf2,
        edge2Cuts
    );

    fileName sFeatFileName
    (
        surf1Name.lessExt().name()
      + "_"
      + surf2Name.lessExt().name()
      + "_"
      + action
    );

    label nFeatEds = inter.cutEdges().size();

    DynamicList<vector> normals(2*nFeatEds);
    vectorField edgeDirections(nFeatEds, vector::zero);
    DynamicList<extendedFeatureEdgeMesh::sideVolumeType> normalVolumeTypes
    (
        2*nFeatEds
    );
    List<DynamicList<label> > edgeNormals(nFeatEds);
    List<DynamicList<label> > normalDirections(nFeatEds);

    forAllConstIter(labelPairLookup, inter.facePairToEdge(), iter)
    {
        const label& cutEdgeI = iter();
        const labelPair& facePair = iter.key();

        const edge& fE = inter.cutEdges()[cutEdgeI];

        const vector& norm1 = surf1.faceNormals()[facePair.first()];
        const vector& norm2 = surf2.faceNormals()[facePair.second()];

        DynamicList<label>& eNormals = edgeNormals[cutEdgeI];
        DynamicList<label>& nDirections = normalDirections[cutEdgeI];

        edgeDirections[cutEdgeI] = fE.vec(inter.cutPoints());

        normals.append(norm1);
        eNormals.append(normals.size() - 1);

        if (surf1Baffle)
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

            nDirections.append(1);
        }
        else
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);
            nDirections.append
            (
                calcNormalDirection
                (
                    norm1,
                    norm2,
                    edgeDirections[cutEdgeI],
                    surf1[facePair.first()].centre(surf1.points()),
                    inter.cutPoints()[fE.start()]
                )
            );
        }

        normals.append(norm2);
        eNormals.append(normals.size() - 1);

        if (surf2Baffle)
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

            nDirections.append(1);
        }
        else
        {
            normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);

            nDirections.append
            (
                calcNormalDirection
                (
                    norm2,
                    norm1,
                    edgeDirections[cutEdgeI],
                    surf2[facePair.second()].centre(surf2.points()),
                    inter.cutPoints()[fE.start()]
                )
            );
        }


        if (surf1Baffle)
        {
            normals.append(norm2);

            if (surf2Baffle)
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

                nDirections.append(1);
            }
            else
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);

                nDirections.append
                (
                    calcNormalDirection
                    (
                        norm2,
                        norm1,
                        edgeDirections[cutEdgeI],
                        surf2[facePair.second()].centre(surf2.points()),
                        inter.cutPoints()[fE.start()]
                    )
                );
            }

            eNormals.append(normals.size() - 1);
        }

        if (surf2Baffle)
        {
            normals.append(norm1);

            if (surf1Baffle)
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::BOTH);

                nDirections.append(1);
            }
            else
            {
                normalVolumeTypes.append(extendedFeatureEdgeMesh::INSIDE);

                nDirections.append
                (
                    calcNormalDirection
                    (
                        norm1,
                        norm2,
                        edgeDirections[cutEdgeI],
                        surf1[facePair.first()].centre(surf1.points()),
                        inter.cutPoints()[fE.start()]
                    )
                );
            }

            eNormals.append(normals.size() - 1);
        }
    }


    label internalStart = -1;
    label nIntOrExt = 0;
    label nFlat = 0;
    label nOpen = 0;
    label nMultiple = 0;

    forAll(edgeNormals, eI)
    {
        label nEdNorms = edgeNormals[eI].size();

        if (nEdNorms == 1)
        {
            nOpen++;
        }
        else if (nEdNorms == 2)
        {
            const vector& n0(normals[edgeNormals[eI][0]]);
            const vector& n1(normals[edgeNormals[eI][1]]);

            if ((n0 & n1) > extendedFeatureEdgeMesh::cosNormalAngleTol_)
            {
                nFlat++;
            }
            else
            {
                nIntOrExt++;
            }
        }
        else if (nEdNorms > 2)
        {
            nMultiple++;
        }
    }

    if (validActions[action] == booleanSurface::UNION)
    {
        if (!invertedSpace)
        {
            // All edges are internal
            internalStart = 0;
        }
        else
        {
            // All edges are external
            internalStart = nIntOrExt;
        }
    }
    else if (validActions[action] == booleanSurface::INTERSECTION)
    {
        if (!invertedSpace)
        {
            // All edges are external
            internalStart = nIntOrExt;
        }
        else
        {
            // All edges are internal
            internalStart = 0;
        }
    }
    else if (validActions[action] == booleanSurface::DIFFERENCE)
    {
        // All edges are external
        internalStart = nIntOrExt;
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Unsupported booleanSurface:booleanOpType and space "
            << action << " " << invertedSpace
            << abort(FatalError);
    }

    // There are no feature points supported by surfaceIntersection
    // Flat, open or multiple edges are assumed to be impossible
    // Region edges are not explicitly supported by surfaceIntersection

    vectorField normalsTmp(normals);
    List<extendedFeatureEdgeMesh::sideVolumeType> normalVolumeTypesTmp
    (
        normalVolumeTypes
    );
    labelListList edgeNormalsTmp(edgeNormals.size());
    forAll(edgeNormalsTmp, i)
    {
        edgeNormalsTmp[i] = edgeNormals[i];
    }
    labelListList normalDirectionsTmp(normalDirections.size());
    forAll(normalDirectionsTmp, i)
    {
        normalDirectionsTmp[i] = normalDirections[i];
    }

    calcFeaturePoints(inter.cutPoints(), inter.cutEdges());

    extendedFeatureEdgeMesh feMesh
    (
        IOobject
        (
            sFeatFileName + ".extendedFeatureEdgeMesh",
            runTime.constant(),
            "extendedFeatureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inter.cutPoints(),
        inter.cutEdges(),

        0,                  // concaveStart,
        0,                  // mixedStart,
        0,                  // nonFeatureStart,

        internalStart,      // internalStart,
        nIntOrExt,           // flatStart,
        nIntOrExt + nFlat,   // openStart,
        nIntOrExt + nFlat + nOpen,   // multipleStart,

        normalsTmp,
        normalVolumeTypesTmp,
        edgeDirections,
        normalDirectionsTmp,
        edgeNormalsTmp,

        labelListList(0),   // featurePointNormals,
        labelListList(0),   // featurePointEdges,
        labelList(0)        // regionEdges
    );

    feMesh.write();

    feMesh.writeObj(feMesh.path()/sFeatFileName);

    {
        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                sFeatFileName + ".eMesh",   // name
                runTime.constant(),                         // instance
                "triSurface",
                runTime,                                    // registry
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.objectPath() << endl;

        bfeMesh.regIOobject::write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
