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
    syncToolsTest

Description
    Test some functionality in syncTools.

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "Random.H"
#include "PackedList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void testPackedList(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing PackedList synchronisation." << endl;

    {
        PackedList<3> bits(mesh.nEdges());
        forAll(bits, i)
        {
            bits.set(i, rndGen.sampleAB<label>(0, 4));
        }

        labelList edgeValues(mesh.nEdges());
        forAll(bits, i)
        {
            edgeValues[i] = bits.get(i);
        }

        PackedList<3> maxBits(bits);
        labelList maxEdgeValues(edgeValues);

        syncTools::syncEdgeList(mesh, bits, minEqOp<unsigned int>(), 0);
        syncTools::syncEdgeList(mesh, edgeValues, minEqOp<label>(), 0);

        syncTools::syncEdgeList(mesh, maxBits, maxEqOp<unsigned int>(), 0);
        syncTools::syncEdgeList
        (
            mesh,
            maxEdgeValues,
            maxEqOp<label>(),
            0
        );

        forAll(bits, i)
        {
            if
            (
                edgeValues[i] != label(bits.get(i))
             || maxEdgeValues[i] != label(maxBits.get(i))
            )
            {
                FatalErrorInFunction
                    << "edge:" << i
                    << " minlabel:" << edgeValues[i]
                    << " minbits:" << bits.get(i)
                    << " maxLabel:" << maxEdgeValues[i]
                    << " maxBits:" << maxBits.get(i)
                    << exit(FatalError);
            }
        }
    }

    {
        PackedList<3> bits(mesh.nPoints());
        forAll(bits, i)
        {
            bits.set(i, rndGen.sampleAB<label>(0, 4));
        }

        labelList pointValues(mesh.nPoints());
        forAll(bits, i)
        {
            pointValues[i] = bits.get(i);
        }

        PackedList<3> maxBits(bits);
        labelList maxPointValues(pointValues);

        syncTools::syncPointList(mesh, bits, minEqOp<unsigned int>(), 0);
        syncTools::syncPointList(mesh, pointValues, minEqOp<label>(), 0);

        syncTools::syncPointList(mesh, maxBits, maxEqOp<unsigned int>(), 0);
        syncTools::syncPointList
        (
            mesh,
            maxPointValues,
            maxEqOp<label>(),
            0
        );

        forAll(bits, i)
        {
            if
            (
                pointValues[i] != label(bits.get(i))
             || maxPointValues[i] != label(maxBits.get(i))
            )
            {
                FatalErrorInFunction
                    << "point:" << i
                    << " at:" << mesh.points()[i]
                    << " minlabel:" << pointValues[i]
                    << " minbits:" << bits.get(i)
                    << " maxLabel:" << maxPointValues[i]
                    << " maxBits:" << maxBits.get(i)
                    << exit(FatalError);
            }
        }
    }

    {
        PackedList<3> bits(mesh.nFaces());
        forAll(bits, facei)
        {
            bits.set(facei, rndGen.sampleAB<label>(0, 4));
        }

        labelList faceValues(mesh.nFaces());
        forAll(bits, facei)
        {
            faceValues[facei] = bits.get(facei);
        }

        PackedList<3> maxBits(bits);
        labelList maxFaceValues(faceValues);

        syncTools::syncFaceList(mesh, bits, minEqOp<unsigned int>());
        syncTools::syncFaceList(mesh, faceValues, minEqOp<label>());

        syncTools::syncFaceList(mesh, maxBits, maxEqOp<unsigned int>());
        syncTools::syncFaceList(mesh, maxFaceValues, maxEqOp<label>());

        forAll(bits, facei)
        {
            if
            (
                faceValues[facei] != label(bits.get(facei))
             || maxFaceValues[facei] != label(maxBits.get(facei))
            )
            {
                FatalErrorInFunction
                    << "face:" << facei
                    << " minlabel:" << faceValues[facei]
                    << " minbits:" << bits.get(facei)
                    << " maxLabel:" << maxFaceValues[facei]
                    << " maxBits:" << maxBits.get(facei)
                    << exit(FatalError);
            }
        }
    }
}


void testSparseData(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing Map synchronisation." << endl;

    WarningInFunction
        << "Position test of sparse data only correct for cases without cyclics"
        << " with shared points." << endl;

    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nFaces()-mesh.nInternalFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );
    const pointField& localPoints = allBoundary.localPoints();


    // Point data
    // ~~~~~~~~~~

    {
        // Create some data. Use slightly perturbed positions.
        Map<point> sparseData;
        pointField fullData(mesh.nPoints(), point(great, great, great));

        forAll(localPoints, i)
        {
            const point pt = localPoints[i] + 1e-4*rndGen.sample01<vector>();

            label meshPointi = allBoundary.meshPoints()[i];

            sparseData.insert(meshPointi, pt);
            fullData[meshPointi] = pt;
        }

        // Pout<< "sparseData:" << sparseData << endl;

        syncTools::syncPointMap
        (
            mesh,
            sparseData,
            minMagSqrEqOp<point>()
            // true                    // apply separation
        );
        syncTools::syncPointList
        (
            mesh,
            fullData,
            minMagSqrEqOp<point>(),
            point(great, great, great)
            // true                    // apply separation
        );

        // Compare.
        // 1. Is all fullData also present in sparseData and same value
        forAll(fullData, meshPointi)
        {
            const point& fullPt = fullData[meshPointi];

            if (fullPt != point(great, great, great))
            {
                const point& sparsePt = sparseData[meshPointi];

                if (fullPt != sparsePt)
                {
                    FatalErrorInFunction
                        << "point:" << meshPointi
                        << " full:" << fullPt
                        << " sparse:" << sparsePt
                        << exit(FatalError);
                }
            }
        }

        // 2. Does sparseData contain more?
        forAllConstIter(Map<point>, sparseData, iter)
        {
            const point& sparsePt = iter();
            label meshPointi = iter.key();
            const point& fullPt = fullData[meshPointi];

            if (fullPt != sparsePt)
            {
                FatalErrorInFunction
                    << "point:" << meshPointi
                    << " full:" << fullPt
                    << " sparse:" << sparsePt
                    << exit(FatalError);
            }
        }
    }


    // Edge data
    // ~~~~~~~~~

    {
        // Create some data. Use slightly perturbed positions.
        EdgeMap<point> sparseData;
        pointField fullData(mesh.nEdges(), point(great, great, great));

        const edgeList& edges = allBoundary.edges();
        const labelList meshEdges = allBoundary.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        );

        forAll(edges, i)
        {
            const edge& e = edges[i];

            const point pt =
                e.centre(localPoints) + 1e-4*rndGen.sample01<vector>();

            label meshEdgeI = meshEdges[i];

            sparseData.insert(mesh.edges()[meshEdgeI], pt);
            fullData[meshEdgeI] = pt;
        }

        // Pout<< "sparseData:" << sparseData << endl;

        syncTools::syncEdgeMap
        (
            mesh,
            sparseData,
            minMagSqrEqOp<point>()
        );
        syncTools::syncEdgeList
        (
            mesh,
            fullData,
            minMagSqrEqOp<point>(),
            point(great, great, great)
        );

        // Compare.
        // 1. Is all fullData also present in sparseData and same value
        forAll(fullData, meshEdgeI)
        {
            const point& fullPt = fullData[meshEdgeI];

            if (fullPt != point(great, great, great))
            {
                const point& sparsePt = sparseData[mesh.edges()[meshEdgeI]];

                if (fullPt != sparsePt)
                {
                    FatalErrorInFunction
                        << "edge:" << meshEdgeI
                        << " points:" << mesh.edges()[meshEdgeI]
                        << " full:" << fullPt
                        << " sparse:" << sparsePt
                        << exit(FatalError);
                }
            }
        }

        // 2. Does sparseData contain more?
        forAll(fullData, meshEdgeI)
        {
            const edge& e = mesh.edges()[meshEdgeI];

            EdgeMap<point>::const_iterator iter = sparseData.find(e);

            if (iter != sparseData.end())
            {
                const point& sparsePt = iter();
                const point& fullPt = fullData[meshEdgeI];

                if (fullPt != sparsePt)
                {
                    FatalErrorInFunction
                        << "Extra edge:" << meshEdgeI
                        << " points:" << mesh.edges()[meshEdgeI]
                        << " full:" << fullPt
                        << " sparse:" << sparsePt
                        << exit(FatalError);
                }
            }
        }
    }
}


void testPointSync(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing point-wise data synchronisation." << endl;

    // Test position.

    {
        pointField syncedPoints(mesh.points());
        syncTools::syncPointPositions
        (
            mesh,
            syncedPoints,
            minMagSqrEqOp<point>(),
            point(great, great, great)
        );

        forAll(syncedPoints, pointi)
        {
            if (mag(syncedPoints[pointi] - mesh.points()[pointi]) > small)
            {
                FatalErrorInFunction
                    << "Point " << pointi
                    << " original location " << mesh.points()[pointi]
                    << " synced location " << syncedPoints[pointi]
                    << exit(FatalError);
            }
        }
    }

    // Test masterPoints

    {
        labelList nMasters(mesh.nPoints(), 0);

        PackedBoolList isMasterPoint(syncTools::getMasterPoints(mesh));

        forAll(isMasterPoint, pointi)
        {
            if (isMasterPoint[pointi])
            {
                nMasters[pointi] = 1;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nMasters,
            plusEqOp<label>(),
            0
        );

        forAll(nMasters, pointi)
        {
            if (nMasters[pointi] != 1)
            {
                WarningInFunction
                    << "Point " << pointi
                    << " original location " << mesh.points()[pointi]
                    << " has " << nMasters[pointi]
                    << " masters."
                    << endl;
            }
        }
    }
}


void testEdgeSync(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing edge-wise data synchronisation." << endl;

    const edgeList& edges = mesh.edges();

    // Test position.

    {
        pointField syncedMids(edges.size());
        forAll(syncedMids, edgeI)
        {
            syncedMids[edgeI] = edges[edgeI].centre(mesh.points());
        }
        syncTools::syncEdgePositions
        (
            mesh,
            syncedMids,
            minMagSqrEqOp<point>(),
            point(great, great, great)
        );

        forAll(syncedMids, edgeI)
        {
            point eMid = edges[edgeI].centre(mesh.points());

            if (mag(syncedMids[edgeI] - eMid) > small)
            {
                FatalErrorInFunction
                    << "Edge " << edgeI
                    << " original midpoint " << eMid
                    << " synced location " << syncedMids[edgeI]
                    << exit(FatalError);
            }
        }
    }

    // Test masterEdges

    {
        labelList nMasters(edges.size(), 0);

        PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));

        forAll(isMasterEdge, edgeI)
        {
            if (isMasterEdge[edgeI])
            {
                nMasters[edgeI] = 1;
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            nMasters,
            plusEqOp<label>(),
            0
        );

        forAll(nMasters, edgeI)
        {
            if (nMasters[edgeI] != 1)
            {
                const edge& e = edges[edgeI];
                WarningInFunction
                    << "Edge " << edgeI
                    << " at:" << mesh.points()[e[0]] << mesh.points()[e[1]]
                    << " has " << nMasters[edgeI]
                    << " masters."
                    << endl;
            }
        }
    }
}


void testFaceSync(const polyMesh& mesh, Random& rndGen)
{
    Info<< nl << "Testing face-wise data synchronisation." << endl;

    // Test position.

    {
        pointField syncedFc(mesh.faceCentres());

        syncTools::syncFacePositions
        (
            mesh,
            syncedFc,
            maxMagSqrEqOp<point>()
        );

        forAll(syncedFc, facei)
        {
            if (mag(syncedFc[facei] - mesh.faceCentres()[facei]) > small)
            {
                FatalErrorInFunction
                    << "Face " << facei
                    << " original centre " << mesh.faceCentres()[facei]
                    << " synced centre " << syncedFc[facei]
                    << exit(FatalError);
            }
        }
    }

    // Test masterFaces

    {
        labelList nMasters(mesh.nFaces(), 0);

        PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

        forAll(isMasterFace, facei)
        {
            if (isMasterFace[facei])
            {
                nMasters[facei] = 1;
            }
        }

        syncTools::syncFaceList
        (
            mesh,
            nMasters,
            plusEqOp<label>()
        );

        forAll(nMasters, facei)
        {
            if (nMasters[facei] != 1)
            {
                FatalErrorInFunction
                    << "Face " << facei
                    << " centre " << mesh.faceCentres()[facei]
                    << " has " << nMasters[facei]
                    << " masters."
                    << exit(FatalError);
            }
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"


    Random rndGen(5341*(Pstream::myProcNo()+1));


    // Face sync
    testFaceSync(mesh, rndGen);

    // Edge sync
    testEdgeSync(mesh, rndGen);

    // Point sync
    testPointSync(mesh, rndGen);

    // PackedList synchronisation
    testPackedList(mesh, rndGen);

    // Sparse synchronisation
    testSparseData(mesh, rndGen);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
