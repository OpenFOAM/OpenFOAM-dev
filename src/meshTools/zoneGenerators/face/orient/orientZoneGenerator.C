/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "orientZoneGenerator.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "patchFaceOrientation.H"
#include "PatchEdgeFaceWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(orientZoneGenerator, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            orientZoneGenerator,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::boolList Foam::zoneGenerators::orientZoneGenerator::normalOrientation
(
    const labelList& faceIndices
) const
{
    const vectorField& faceAreas = mesh_.faceAreas();

    boolList flipMap(faceIndices.size(), false);

    forAll(faceIndices, fi)
    {
        const label facei = faceIndices[fi];

        if ((faceAreas[facei] & normal_) < 0)
        {
            flipMap[fi] = true;
        }
    }

    return flipMap;
}


Foam::boolList Foam::zoneGenerators::orientZoneGenerator::pointOrientation
(
    const faceZone& fZone
) const
{
    if (fZone.checkParallelSync())
    {
        FatalErrorInFunction
            << "Face zone " << fZone.name()
            << " is not parallel synchronised."
            << " Any coupled face also needs its coupled version to be included"
            << " and with opposite flipMap."
            << exit(FatalError);
    }

    const labelList& faceIndices = fZone;

    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), faceIndices),
        mesh_.points()
    );

    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    // Data on all edges and faces
    List<patchFaceOrientation> allEdgeInfo(patch.nEdges());
    List<patchFaceOrientation> allFaceInfo(patch.size());

    // Make sure we don't walk through
    // - slaves of coupled faces
    // - non-manifold edges
    {
        const polyBoundaryMesh& bm = mesh_.boundaryMesh();

        label nProtected = 0;

        forAll(faceIndices, facei)
        {
            const label meshFacei = faceIndices[facei];
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

        if (debug)
        {
            Info<< "Protected from visiting "
                << returnReduce(nProtected, sumOp<label>())
                << " slaves of coupled faces" << nl << endl;
        }
    }
    {
        // Number of (master)faces per edge
        labelList nMasterFaces(patch.nEdges(), 0);

        forAll(faceIndices, facei)
        {
            const label meshFacei = faceIndices[facei];

            if (isMasterFace[meshFacei])
            {
                const labelList& fEdges = patch.faceEdges()[facei];
                forAll(fEdges, fEdgei)
                {
                    nMasterFaces[fEdges[fEdgei]]++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            patch.meshEdges(mesh_.edges(), mesh_.pointEdges()),
            nMasterFaces,
            plusEqOp<label>(),
            label(0)
        );


        label nProtected = 0;

        forAll(nMasterFaces, edgei)
        {
            if (nMasterFaces[edgei] > 2)
            {
                allEdgeInfo[edgei] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        if (debug)
        {
            Info<< "Protected from visiting "
                << returnReduce(nProtected, sumOp<label>())
                << " non-manifold edges" << nl << endl;
        }
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

        const label proci = globalFaces.whichProcID(unsetFacei);
        const label seedFacei = globalFaces.toLocal(proci, unsetFacei);

        if (debug)
        {
            Info<< "Seeding from processor " << proci << " face " << seedFacei
                << endl;
        }

        if (proci == Pstream::myProcNo())
        {
            // Determine orientation of seedFace

            const vector d = outsidePoint_ - patch.faceCentres()[seedFacei];
            const vector& fn = patch.faceNormals()[seedFacei];

            // Set information to correct orientation
            patchFaceOrientation& faceInfo = allFaceInfo[seedFacei];
            faceInfo = orientedSurface::NOFLIP;

            if ((fn&d) < 0)
            {
                faceInfo.flip();

                if (debug)
                {
                    Pout<< "Face " << seedFacei << " at "
                        << patch.faceCentres()[seedFacei]
                        << " with normal " << fn
                        << " needs to be flipped." << endl;
                }
            }
            else if (debug)
            {
                Pout<< "Face " << seedFacei << " at "
                    << patch.faceCentres()[seedFacei]
                    << " with normal " << fn
                    << " points in positive direction (cos = " << (fn&d)/mag(d)
                    << ")" << endl;
            }


            const labelList& fEdges = patch.faceEdges()[seedFacei];
            forAll(fEdges, fEdgei)
            {
                const label edgei = fEdges[fEdgei];

                patchFaceOrientation& edgeInfo = allEdgeInfo[edgei];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgei,
                        seedFacei,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgei);
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
            mesh_,
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
        const polyBoundaryMesh& bm = mesh_.boundaryMesh();

        labelList neiStatus
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            orientedSurface::UNVISITED
        );

        forAll(faceIndices, i)
        {
            const label meshFacei = faceIndices[i];
            if (!mesh_.isInternalFace(meshFacei))
            {
                neiStatus[meshFacei-mesh_.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiStatus);

        forAll(faceIndices, i)
        {
            const label meshFacei = faceIndices[i];
            const label patchi = bm.whichPatch(meshFacei);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFacei]
            )
            {
                // Slave side. Take flipped from neighbour
                const label bFacei = meshFacei-mesh_.nInternalFaces();

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


    // Convert to flipMap

    boolList flipMap(faceIndices.size());

    forAll(allFaceInfo, facei)
    {
        if (allFaceInfo[facei] == orientedSurface::NOFLIP)
        {
            flipMap[facei] = false;
        }
        else if (allFaceInfo[facei] == orientedSurface::FLIP)
        {
            flipMap[facei] = true;
        }
        else
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << facei
                << " centre:" << mesh_.faceCentres()[faceIndices[facei]]
                << abort(FatalError);
        }
    }

    return flipMap;
}


Foam::boolList Foam::zoneGenerators::orientZoneGenerator::orientation
(
    const faceZone& fZone
) const
{
    if (magSqr(normal_) > 0)
    {
        return normalOrientation(fZone);
    }
    else
    {
        return pointOrientation(fZone);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::orientZoneGenerator::orientZoneGenerator
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneGenerator_(zoneGenerator::New(mesh, dict)),
    normal_(dict.lookupOrDefault<vector>("normal", Zero)),
    outsidePoint_
    (
        !dict.found("normal")
       ? dict.lookup<vector>("outsidePoint")
       : Zero
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::orientZoneGenerator::~orientZoneGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::orientZoneGenerator::generate() const
{
    zoneSet zs(zoneGenerator_->generate());
    const faceZone& fZone = zs.fZone();

    return zoneSet
    (
        new faceZone
        (
            zoneName_,
            fZone,
            orientation(fZone),
            mesh_.faceZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
