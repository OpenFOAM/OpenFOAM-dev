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

Application
    surfaceCheck

Description
    Checks geometric and topological quality of a surface.

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "argList.H"
#include "OFstream.H"
#include "OBJstream.H"
#include "SortableList.H"
#include "PatchTools.H"
#include "vtkSurfaceWriter.H"

using namespace Foam;

// Does face use valid vertices?
bool validTri
(
    const bool verbose,
    const triSurface& surf,
    const label facei
)
{
    // Simple check on indices ok.

    const labelledTri& f = surf[facei];

    forAll(f, fp)
    {
        if (f[fp] < 0 || f[fp] >= surf.points().size())
        {
            WarningInFunction
                << "triangle " << facei << " vertices " << f
                << " uses point indices outside point range 0.."
                << surf.points().size()-1 << endl;
            return false;
        }
    }

    if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
    {
        WarningInFunction
            << "triangle " << facei
            << " uses non-unique vertices " << f
            << " coords:" << f.points(surf.points())
            << endl;
        return false;
    }

    // duplicate triangle check

    const labelList& fFaces = surf.faceFaces()[facei];

    // Check if faceNeighbours use same points as this face.
    // Note: discards normal information - sides of baffle are merged.
    forAll(fFaces, i)
    {
        label nbrFacei = fFaces[i];

        if (nbrFacei <= facei)
        {
            // lower numbered faces already checked
            continue;
        }

        const labelledTri& nbrF = surf[nbrFacei];

        if
        (
            ((f[0] == nbrF[0]) || (f[0] == nbrF[1]) || (f[0] == nbrF[2]))
         && ((f[1] == nbrF[0]) || (f[1] == nbrF[1]) || (f[1] == nbrF[2]))
         && ((f[2] == nbrF[0]) || (f[2] == nbrF[1]) || (f[2] == nbrF[2]))
        )
        {
            WarningInFunction
                << "triangle " << facei << " vertices " << f
                << " has the same vertices as triangle " << nbrFacei
                << " vertices " << nbrF
                << " coords:" << f.points(surf.points())
                << endl;

            return false;
        }
    }
    return true;
}


labelList countBins
(
    const scalar min,
    const scalar max,
    const label nBins,
    const scalarField& vals
)
{
    scalar dist = nBins/(max - min);

    labelList binCount(nBins, 0);

    forAll(vals, i)
    {
        scalar val = vals[i];

        label index = -1;

        if (Foam::mag(val - min) < small)
        {
            index = 0;
        }
        else if (val >= max - small)
        {
            index = nBins - 1;
        }
        else
        {
            index = label((val - min)*dist);

            if ((index < 0) || (index >= nBins))
            {
                WarningInFunction
                    << "value " << val << " at index " << i
                    << " outside range " << min << " .. " << max << endl;

                if (index < 0)
                {
                    index = 0;
                }
                else
                {
                    index = nBins - 1;
                }
            }
        }
        binCount[index]++;
    }

    return binCount;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validArgs.append("surface file");
    argList::addBoolOption
    (
        "checkSelfIntersection",
        "also check for self-intersection"
    );
    argList::addBoolOption
    (
        "splitNonManifold",
        "split surface along non-manifold edges"
        " (default split is fully disconnected)"
    );
    argList::addBoolOption
    (
        "verbose",
        "verbose operation"
    );
    argList::addBoolOption
    (
        "blockMesh",
        "write vertices/blocks for blockMeshDict"
    );

    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const bool checkSelfIntersect = args.optionFound("checkSelfIntersection");
    const bool verbose = args.optionFound("verbose");
    const bool splitNonManifold = args.optionFound("splitNonManifold");

    Info<< "Reading surface from " << surfFileName << " ..." << nl << endl;


    // Read
    // ~~~~

    triSurface surf(surfFileName);


    Info<< "Statistics:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    // write bounding box corners
    if (args.optionFound("blockMesh"))
    {
        pointField cornerPts(boundBox(surf.points(), false).points());

        Info<< "// blockMeshDict info" << nl << nl;

        Info<< "vertices\n(" << nl;
        forAll(cornerPts, ptI)
        {
            Info<< "    " << cornerPts[ptI] << nl;
        }

        // number of divisions needs adjustment later
        Info<< ");\n" << nl
            << "blocks\n"
            << "(\n"
            << "    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)\n"
            << ");\n" << nl;

        Info<< "edges\n();" << nl
            << "patches\n();" << endl;

        Info<< nl << "// end blockMeshDict info" << nl << endl;
    }


    // Region sizes
    // ~~~~~~~~~~~~

    {
        labelList regionSize(surf.patches().size(), 0);

        forAll(surf, facei)
        {
            label region = surf[facei].region();

            if (region < 0 || region >= regionSize.size())
            {
                WarningInFunction
                    << "Triangle " << facei << " vertices " << surf[facei]
                    << " has region " << region << " which is outside the range"
                    << " of regions 0.." << surf.patches().size()-1
                    << endl;
            }
            else
            {
                regionSize[region]++;
            }
        }

        Info<< "Region\tSize" << nl
            << "------\t----" << nl;
        forAll(surf.patches(), patchi)
        {
            Info<< surf.patches()[patchi].name() << '\t'
                << regionSize[patchi] << nl;
        }
        Info<< nl << endl;
    }


    // Check triangles
    // ~~~~~~~~~~~~~~~

    {
        DynamicList<label> illegalFaces(surf.size()/100 + 1);

        forAll(surf, facei)
        {
            if (!validTri(verbose, surf, facei))
            {
                illegalFaces.append(facei);
            }
        }

        if (illegalFaces.size())
        {
            Info<< "Surface has " << illegalFaces.size()
                << " illegal triangles." << endl;

            OFstream str("illegalFaces");
            Info<< "Dumping conflicting face labels to " << str.name() << endl
                << "Paste this into the input for surfaceSubset" << endl;
            str << illegalFaces;
        }
        else
        {
            Info<< "Surface has no illegal triangles." << endl;
        }
        Info<< endl;
    }



    // Triangle quality
    // ~~~~~~~~~~~~~~~~

    {
        scalarField triQ(surf.size(), 0);
        forAll(surf, facei)
        {
            const labelledTri& f = surf[facei];

            if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
            {
                // WarningIn(args.executable())
                //    << "Illegal triangle " << facei << " vertices " << f
                //    << " coords " << f.points(surf.points()) << endl;
            }
            else
            {
                triQ[facei] = triPointRef
                (
                    surf.points()[f[0]],
                    surf.points()[f[1]],
                    surf.points()[f[2]]
                ).quality();
            }
        }

        labelList binCount = countBins(0, 1, 20, triQ);

        Info<< "Triangle quality (equilateral=1, collapsed=0):"
            << endl;


        OSstream& os = Info;
        os.width(4);

        scalar dist = (1.0 - 0.0)/20.0;
        scalar min = 0;
        forAll(binCount, binI)
        {
            Info<< "    " << min << " .. " << min+dist << "  : "
                << 1.0/surf.size() * binCount[binI]
                << endl;
            min += dist;
        }
        Info<< endl;

        label minIndex = findMin(triQ);
        label maxIndex = findMax(triQ);

        Info<< "    min " << triQ[minIndex] << " for triangle " << minIndex
            << nl
            << "    max " << triQ[maxIndex] << " for triangle " << maxIndex
            << nl
            << endl;


        if (triQ[minIndex] < small)
        {
            WarningInFunction
                << triQ[minIndex] << ". This might give problems in"
                << " self-intersection testing later on." << endl;
        }

        // Dump for subsetting
        {
            DynamicList<label> problemFaces(surf.size()/100+1);

            forAll(triQ, facei)
            {
                if (triQ[facei] < 1e-11)
                {
                    problemFaces.append(facei);
                }
            }

            if (!problemFaces.empty())
            {
                OFstream str("badFaces");

                Info<< "Dumping bad quality faces to " << str.name() << endl
                    << "Paste this into the input for surfaceSubset" << nl
                    << nl << endl;

                str << problemFaces;
            }
        }
    }



    // Edges
    // ~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        scalarField edgeMag(edges.size());

        forAll(edges, edgeI)
        {
            edgeMag[edgeI] = edges[edgeI].mag(localPoints);
        }

        label minEdgeI = findMin(edgeMag);
        label maxEdgeI = findMax(edgeMag);

        const edge& minE = edges[minEdgeI];
        const edge& maxE = edges[maxEdgeI];


        Info<< "Edges:" << nl
            << "    min " << edgeMag[minEdgeI] << " for edge " << minEdgeI
            << " points " << localPoints[minE[0]] << localPoints[minE[1]]
            << nl
            << "    max " << edgeMag[maxEdgeI] << " for edge " << maxEdgeI
            << " points " << localPoints[maxE[0]] << localPoints[maxE[1]]
            << nl
            << endl;
    }



    // Close points
    // ~~~~~~~~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        const boundBox bb(localPoints);
        scalar smallDim = 1e-6 * bb.mag();

        Info<< "Checking for points less than 1e-6 of bounding box ("
            << bb.span() << " metre) apart."
            << endl;

        // Sort points
        SortableList<scalar> sortedMag(mag(localPoints));

        label nClose = 0;

        for (label i = 1; i < sortedMag.size(); i++)
        {
            label ptI = sortedMag.indices()[i];

            label prevPtI = sortedMag.indices()[i-1];

            if (mag(localPoints[ptI] - localPoints[prevPtI]) < smallDim)
            {
                // Check if neighbours.
                const labelList& pEdges = surf.pointEdges()[ptI];

                label edgeI = -1;

                forAll(pEdges, i)
                {
                    const edge& e = edges[pEdges[i]];

                    if (e[0] == prevPtI || e[1] == prevPtI)
                    {
                        // point1 and point0 are connected through edge.
                        edgeI = pEdges[i];

                        break;
                    }
                }

                nClose++;

                if (edgeI == -1)
                {
                    Info<< "    close unconnected points "
                        << ptI << ' ' << localPoints[ptI]
                        << " and " << prevPtI << ' '
                        << localPoints[prevPtI]
                        << " distance:"
                        << mag(localPoints[ptI] - localPoints[prevPtI])
                        << endl;
                }
                else
                {
                    Info<< "    small edge between points "
                        << ptI << ' ' << localPoints[ptI]
                        << " and " << prevPtI << ' '
                        << localPoints[prevPtI]
                        << " distance:"
                        << mag(localPoints[ptI] - localPoints[prevPtI])
                        << endl;
                }
            }
        }

        Info<< "Found " << nClose << " nearby points." << nl
            << endl;
    }



    // Check manifold
    // ~~~~~~~~~~~~~~

    DynamicList<label> problemFaces(surf.size()/100 + 1);

    const labelListList& eFaces = surf.edgeFaces();

    label nSingleEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() == 1)
        {
            problemFaces.append(myFaces[0]);

            nSingleEdges++;
        }
    }

    label nMultEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() > 2)
        {
            forAll(myFaces, myFacei)
            {
                problemFaces.append(myFaces[myFacei]);
            }

            nMultEdges++;
        }
    }
    problemFaces.shrink();

    if ((nSingleEdges != 0) || (nMultEdges != 0))
    {
        Info<< "Surface is not closed since not all edges connected to "
            << "two faces:" << endl
            << "    connected to one face : " << nSingleEdges << endl
            << "    connected to >2 faces : " << nMultEdges << endl;

        Info<< "Conflicting face labels:" << problemFaces.size() << endl;

        OFstream str("problemFaces");

        Info<< "Dumping conflicting face labels to " << str.name() << endl
            << "Paste this into the input for surfaceSubset" << endl;

        str << problemFaces;
    }
    else
    {
        Info<< "Surface is closed. All edges connected to two faces." << endl;
    }
    Info<< endl;



    // Check singly connected domain
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        boolList borderEdge(surf.nEdges(), false);
        if (splitNonManifold)
        {
            const labelListList& eFaces = surf.edgeFaces();
            forAll(eFaces, edgeI)
            {
                if (eFaces[edgeI].size() > 2)
                {
                    borderEdge[edgeI] = true;
                }
            }
        }

        labelList faceZone;
        label numZones = surf.markZones(borderEdge, faceZone);

        Info<< "Number of unconnected parts : " << numZones << endl;

        if (numZones > 1)
        {
            Info<< "Splitting surface into parts ..." << endl << endl;

            fileName surfFileNameBase(surfFileName.name());
            const word fileType = surfFileNameBase.ext();
            // Strip extension
            surfFileNameBase = surfFileNameBase.lessExt();
            // If extension was .gz strip original extension
            if (fileType == "gz")
            {
                surfFileNameBase = surfFileNameBase.lessExt();
            }


            {
                Info<< "Writing zoning to "
                    <<  fileName
                        (
                            "zone_"
                          + surfFileNameBase
                          + '.'
                          + vtkSurfaceWriter::typeName
                        )
                    << "..." << endl << endl;

                // Convert data
                scalarField scalarFaceZone(faceZone.size());
                forAll(faceZone, i)
                {
                    scalarFaceZone[i] = faceZone[i];
                }
                faceList faces(surf.size());
                forAll(surf, i)
                {
                    faces[i] = surf[i].triFaceFace();
                }

                vtkSurfaceWriter().write
                (
                    surfFileName.path(),
                    surfFileNameBase,
                    surf.points(),
                    faces,
                    "zone",
                    scalarFaceZone,
                    false               // face based data
                );
            }


            for (label zone = 0; zone < numZones; zone++)
            {
                boolList includeMap(surf.size(), false);

                forAll(faceZone, facei)
                {
                    if (faceZone[facei] == zone)
                    {
                        includeMap[facei] = true;
                    }
                }

                labelList pointMap;
                labelList faceMap;

                triSurface subSurf
                (
                    surf.subsetMesh
                    (
                        includeMap,
                        pointMap,
                        faceMap
                    )
                );

                fileName subName(surfFileNameBase + "_" + name(zone) + ".obj");

                Info<< "writing part " << zone << " size " << subSurf.size()
                    << " to " << subName << endl;

                subSurf.write(subName);
            }
        }
    }



    // Check orientation
    // ~~~~~~~~~~~~~~~~~

    labelHashSet borderEdge(surf.size()/1000);
    PatchTools::checkOrientation(surf, false, &borderEdge);

    //
    // Colour all faces into zones using borderEdge
    //
    labelList normalZone;
    label numNormalZones = PatchTools::markZones(surf, borderEdge, normalZone);

    Info<< endl
        << "Number of zones (connected area with consistent normal) : "
        << numNormalZones << endl;

    if (numNormalZones > 1)
    {
        Info<< "More than one normal orientation." << endl;
    }
    Info<< endl;



    // Check self-intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (checkSelfIntersect)
    {
        Info<< "Checking self-intersection." << endl;

        triSurfaceSearch querySurf(surf);

        const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

        OBJstream intStream("selfInterPoints.obj");

        label nInt = 0;

        forAll(surf.edges(), edgeI)
        {
            const edge& e = surf.edges()[edgeI];

            pointIndexHit hitInfo
            (
                tree.findLine
                (
                    surf.points()[surf.meshPoints()[e[0]]],
                    surf.points()[surf.meshPoints()[e[1]]],
                    treeDataTriSurface::findSelfIntersectOp
                    (
                        tree,
                        edgeI
                    )
                )
            );

            if (hitInfo.hit())
            {
                intStream.write(hitInfo.hitPoint());
                nInt++;
            }
        }

        if (nInt == 0)
        {
            Info<< "Surface is not self-intersecting" << endl;
        }
        else
        {
            Info<< "Surface is self-intersecting at " << nInt
                << " locations." << endl;
            Info<< "Writing intersection points to " << intStream.name()
                << endl;
        }

        // surfaceIntersection inter(querySurf);
        //
        // if (inter.cutEdges().empty() && inter.cutPoints().empty())
        //{
        //    Info<< "Surface is not self-intersecting" << endl;
        //}
        // else
        //{
        //    Info<< "Surface is self-intersecting" << endl;
        //    Info<< "Writing edges of intersection to selfInter.obj" << endl;
        //
        //    OFstream intStream("selfInter.obj");
        //    forAll(inter.cutPoints(), cutPointi)
        //    {
        //        const point& pt = inter.cutPoints()[cutPointi];
        //
        //        intStream << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z()
        //            << endl;
        //    }
        //    forAll(inter.cutEdges(), cutEdgeI)
        //    {
        //        const edge& e = inter.cutEdges()[cutEdgeI];
        //
        //        intStream << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
        //    }
        //}
        Info<< endl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
