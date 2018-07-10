/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    orientFaceZone

Description
    Corrects orientation of faceZone.

    - correct in parallel - excludes coupled faceZones from walk
    - correct for non-manifold faceZones - restarts walk

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "syncTools.H"
#include "patchFaceOrientation.H"
#include "PatchEdgeFaceWave.H"
#include "orientedSurface.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    argList::validArgs.append("faceZone");
    argList::validArgs.append("outsidePoint");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedPolyMesh.H"

    const word zoneName  = args[1];
    const point outsidePoint = args.argRead<point>(2);

    Info<< "Orienting faceZone " << zoneName
        << " such that " << outsidePoint << " is outside"
        << nl << endl;


    const faceZone& fZone = mesh.faceZones()[zoneName];

    if (fZone.checkParallelSync())
    {
        FatalErrorInFunction
            << "Face zone " << fZone.name()
            << " is not parallel synchronised."
            << " Any coupled face also needs its coupled version to be included"
            << " and with opposite flipMap."
            << exit(FatalError);
    }

    const labelList& faceLabels = fZone;

    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );



    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));


    // Data on all edges and faces
    List<patchFaceOrientation> allEdgeInfo(patch.nEdges());
    List<patchFaceOrientation> allFaceInfo(patch.size());

    // Make sure we don't walk through
    // - slaves of coupled faces
    // - non-manifold edges
    {
        const polyBoundaryMesh& bm = mesh.boundaryMesh();

        label nProtected = 0;

        forAll(faceLabels, facei)
        {
            const label meshFacei = faceLabels[facei];
            const label patchi = bm.whichPatch(meshFacei);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFacei]
            )
            {
                // Slave side. Mark so doesn't get visited.
                allFaceInfo[facei] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " slaves of coupled faces" << nl << endl;
    }
    {
        // Number of (master)faces per edge
        labelList nMasterFaces(patch.nEdges(), 0);

        forAll(faceLabels, facei)
        {
            const label meshFacei = faceLabels[facei];

            if (isMasterFace[meshFacei])
            {
                const labelList& fEdges = patch.faceEdges()[facei];
                forAll(fEdges, fEdgeI)
                {
                    nMasterFaces[fEdges[fEdgeI]]++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            patch.meshEdges(mesh.edges(), mesh.pointEdges()),
            nMasterFaces,
            plusEqOp<label>(),
            label(0)
        );


        label nProtected = 0;

        forAll(nMasterFaces, edgeI)
        {
            if (nMasterFaces[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " non-manifold edges" << nl << endl;
    }



    DynamicList<label> changedEdges;
    DynamicList<patchFaceOrientation> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchFaceOrientation
    >::propagationTol();

    int dummyTrackData;

    globalIndex globalFaces(patch.size());

    while (true)
    {
        // Pick an unset face
        label unsetFacei = labelMax;
        forAll(allFaceInfo, facei)
        {
            if (allFaceInfo[facei] == orientedSurface::UNVISITED)
            {
                unsetFacei = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(unsetFacei, minOp<label>());

        if (unsetFacei == labelMax)
        {
            break;
        }

        label proci = globalFaces.whichProcID(unsetFacei);
        label seedFacei = globalFaces.toLocal(proci, unsetFacei);
        Info<< "Seeding from processor " << proci << " face " << seedFacei
            << endl;

        if (proci == Pstream::myProcNo())
        {
            // Determine orientation of seedFace

            vector d = outsidePoint-patch.faceCentres()[seedFacei];
            const vector& fn = patch.faceNormals()[seedFacei];

            // Set information to correct orientation
            patchFaceOrientation& faceInfo = allFaceInfo[seedFacei];
            faceInfo = orientedSurface::NOFLIP;

            if ((fn&d) < 0)
            {
                faceInfo.flip();

                Pout<< "Face " << seedFacei << " at "
                    << patch.faceCentres()[seedFacei]
                    << " with normal " << fn
                    << " needs to be flipped." << endl;
            }
            else
            {
                Pout<< "Face " << seedFacei << " at "
                    << patch.faceCentres()[seedFacei]
                    << " with normal " << fn
                    << " points in positive direction (cos = " << (fn&d)/mag(d)
                    << ")" << endl;
            }


            const labelList& fEdges = patch.faceEdges()[seedFacei];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchFaceOrientation& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh,
                        patch,
                        edgeI,
                        seedFacei,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduce(changedEdges.size(), sumOp<label>()) == 0)
        {
            break;
        }



        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchFaceOrientation
        > calc
        (
            mesh,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );
    }


    // Push master zone info over to slave (since slave faces never visited)
    {
        const polyBoundaryMesh& bm = mesh.boundaryMesh();

        labelList neiStatus
        (
            mesh.nFaces()-mesh.nInternalFaces(),
            orientedSurface::UNVISITED
        );

        forAll(faceLabels, i)
        {
            const label meshFacei = faceLabels[i];
            if (!mesh.isInternalFace(meshFacei))
            {
                neiStatus[meshFacei-mesh.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh, neiStatus);

        forAll(faceLabels, i)
        {
            const label meshFacei = faceLabels[i];
            const label patchi = bm.whichPatch(meshFacei);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFacei]
            )
            {
                // Slave side. Take flipped from neighbour
                label bFacei = meshFacei-mesh.nInternalFaces();

                if (neiStatus[bFacei] == orientedSurface::NOFLIP)
                {
                    allFaceInfo[i] = orientedSurface::FLIP;
                }
                else if (neiStatus[bFacei] == orientedSurface::FLIP)
                {
                    allFaceInfo[i] = orientedSurface::NOFLIP;
                }
                else
                {
                    FatalErrorInFunction
                        << "Incorrect status for face " << meshFacei
                        << abort(FatalError);
                }
            }
        }
    }


    // Convert to flipmap and adapt faceZones

    boolList newFlipMap(allFaceInfo.size(), false);
    label nChanged = 0;
    forAll(allFaceInfo, facei)
    {
        if (allFaceInfo[facei] == orientedSurface::NOFLIP)
        {
            newFlipMap[facei] = false;
        }
        else if (allFaceInfo[facei] == orientedSurface::FLIP)
        {
            newFlipMap[facei] = true;
        }
        else
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << facei
                << " centre:" << mesh.faceCentres()[faceLabels[facei]]
                << abort(FatalError);
        }

        if (fZone.flipMap()[facei] != newFlipMap[facei])
        {
            nChanged++;
        }
    }

    reduce(nChanged, sumOp<label>());
    if (nChanged > 0)
    {
        Info<< "Flipping " << nChanged << " out of "
            << globalFaces.size() << " faces." << nl << endl;

        mesh.faceZones()[zoneName].resetAddressing(faceLabels, newFlipMap);
        if (!mesh.faceZones().write())
        {
            FatalErrorInFunction
                << "Failed writing faceZones" << exit(FatalError);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
