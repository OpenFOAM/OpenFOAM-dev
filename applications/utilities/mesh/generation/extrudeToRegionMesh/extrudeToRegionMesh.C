/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    extrudeToRegionMesh

Description
    Extrude faceZones (internal or boundary faces) or faceSets (boundary faces
    only) into a separate mesh (as a different region).

    - used to e.g. extrude baffles (extrude internal faces) or create
    liquid film regions.
    - if extruding internal faces:
        - create baffles in original mesh with mappedWall patches
    - if extruding boundary faces:
        - convert boundary faces to mappedWall patches
    - extrude edges of faceZone as a \<zone\>_sidePatch
    - extrudes into master direction (i.e. away from the owner cell
      if flipMap is false)

\verbatim

Internal face extrusion
-----------------------

    +-------------+
    |             |
    |             |
    +---AAAAAAA---+
    |             |
    |             |
    +-------------+

    AAA=faceZone to extrude.


For the case of no flipMap the extrusion starts at owner and extrudes
into the space of the neighbour:

      +CCCCCCC+
      |       |         <= extruded mesh
      +BBBBBBB+

    +-------------+
    |             |
    | (neighbour) |
    |___CCCCCCC___|       <= original mesh (with 'baffles' added)
    |   BBBBBBB   |
    |(owner side) |
    |             |
    +-------------+

    BBB=mapped between owner on original mesh and new extrusion.
        (zero offset)
    CCC=mapped between neighbour on original mesh and new extrusion
        (offset due to the thickness of the extruded mesh)

For the case of flipMap the extrusion is the other way around: from the
neighbour side into the owner side.


Boundary face extrusion
-----------------------

    +--AAAAAAA--+
    |           |
    |           |
    +-----------+

    AAA=faceZone to extrude. E.g. slave side is owner side (no flipmap)

becomes

      +CCCCCCC+
      |       |         <= extruded mesh
      +BBBBBBB+

    +--BBBBBBB--+
    |           |       <= original mesh
    |           |
    +-----------+

    BBB=mapped between original mesh and new extrusion
    CCC=polypatch


    Notes:
      - when extruding cyclics with only one cell in between it does not
        detect this as a cyclic since the face is the same face. It will
        only work if the coupled edge extrudes a different face so if there
        are more than 1 cell in between.

\endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyTopoChange.H"
#include "mappedWallPolyPatch.H"
#include "mappedExtrudedWallPolyPatch.H"
#include "createShellMesh.H"
#include "syncTools.H"
#include "cyclicPolyPatch.H"
#include "wedgePolyPatch.H"
#include "extrudeModel.H"
#include "faceSet.H"
#include "fvMeshTools.H"
#include "OBJstream.H"
#include "PatchTools.H"
#include "systemDict.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label findPatchID(const List<polyPatch*>& newPatches, const word& name)
{
    forAll(newPatches, i)
    {
        if (newPatches[i]->name() == name)
        {
            return i;
        }
    }
    return -1;
}


label addPatch
(
    const polyBoundaryMesh& patches,
    const word& patchName,
    const word& patchType,
    const dictionary& dict,
    DynamicList<polyPatch*>& newPatches
)
{
    label patchi = findPatchID(newPatches, patchName);

    if (patchi != -1)
    {
        if (newPatches[patchi]->type() == patchType)
        {
            return patchi;
        }
        else
        {
            FatalErrorInFunction
                << "Already have patch " << patchName
                << " but of type " << newPatches[patchi]->type()
                << exit(FatalError);
        }
    }

    patchi = newPatches.size();

    label startFacei = 0;
    if (patchi > 0)
    {
        const polyPatch& pp = *newPatches.last();
        startFacei = pp.start()+pp.size();
    }

    dictionary patchDict(dict);
    patchDict.set("type", patchType);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", startFacei);

    newPatches.append
    (
        polyPatch::New
        (
            patchName,
            patchDict,
            patchi,
            patches
        ).ptr()
    );

    return patchi;
}


// Remove zero-sized patches
void deleteEmptyPatches(fvMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    wordList masterNames;
    if (Pstream::master())
    {
        masterNames = patches.names();
    }
    Pstream::scatter(masterNames);


    labelList oldToNew(patches.size(), -1);
    label usedI = 0;
    label notUsedI = patches.size();

    // Add all the non-empty, non-processor patches
    forAll(masterNames, masterI)
    {
        label patchi = patches.findPatchID(masterNames[masterI]);

        if (patchi != -1)
        {
            if (isA<processorPolyPatch>(patches[patchi]))
            {
                // Similar named processor patch? Not 'possible'.
                if (patches[patchi].size() == 0)
                {
                    Pout<< "Deleting processor patch " << patchi
                        << " name:" << patches[patchi].name()
                        << endl;
                    oldToNew[patchi] = --notUsedI;
                }
                else
                {
                    oldToNew[patchi] = usedI++;
                }
            }
            else
            {
                // Common patch.
                if (returnReduce(patches[patchi].size(), sumOp<label>()) == 0)
                {
                    Pout<< "Deleting patch " << patchi
                        << " name:" << patches[patchi].name()
                        << endl;
                    oldToNew[patchi] = --notUsedI;
                }
                else
                {
                    oldToNew[patchi] = usedI++;
                }
            }
        }
    }

    // Add remaining patches at the end
    forAll(patches, patchi)
    {
        if (oldToNew[patchi] == -1)
        {
            // Unique to this processor. Note: could check that these are
            // only processor patches.
            if (patches[patchi].size() == 0)
            {
                Pout<< "Deleting processor patch " << patchi
                    << " name:" << patches[patchi].name()
                    << endl;
                oldToNew[patchi] = --notUsedI;
            }
            else
            {
                oldToNew[patchi] = usedI++;
            }
        }
    }

    fvMeshTools::reorderPatches(mesh, oldToNew, usedI, true);
}


enum class connectivity
{
    unset,
    boundary,
    internal,
    error
};


struct isInternalEqOp
{
    void operator()(connectivity& x, const connectivity& y) const
    {
        if (x == connectivity::unset)
        {
            // Set x if not yet set
            x = y;
        }
        else if (y == connectivity::unset)
        {
            // Do not set x if y is not yet set
        }
        else if (x == y)
        {
            // Do nothing if x and y are equal
        }
        else
        {
            // Set an error value if x and y are both set and unequal
            x = connectivity::error;
        }
    }
};


Ostream& operator<<(Ostream& os, const connectivity& p)
{
    return os << static_cast<label>(p);
}


Istream& operator>>(Istream& is, connectivity& p)
{
    label l;
    is >> l;
    p = static_cast<connectivity>(l);
    return is;
}


// Return a bool list indicating whether or not the zone contains internal or
// boundary faces
boolList calcZoneIsInternal
(
    const polyMesh& mesh,
    const wordList& zoneNames,
    const labelList& extrudeFaces,
    const labelList& extrudeFaceZoneIDs
)
{
    List<connectivity> zoneConnectivity(zoneNames.size(), connectivity::unset);

    forAll(extrudeFaces, i)
    {
        const label facei = extrudeFaces[i];
        const label zonei = extrudeFaceZoneIDs[i];
        isInternalEqOp()
        (
            zoneConnectivity[zonei],
            mesh.isInternalFace(facei)
          ? connectivity::internal
          : connectivity::boundary
        );
    }

    Pstream::listCombineGather(zoneConnectivity, isInternalEqOp());
    Pstream::listCombineScatter(zoneConnectivity);

    boolList zoneIsInternal(zoneNames.size());
    forAll(zoneConnectivity, zonei)
    {
        switch (zoneConnectivity[zonei])
        {
            case connectivity::unset:
                // This zone doesn't contain any faces, so it doesn't matter
                // what we do here
                zoneIsInternal[zonei] = false;
                break;
            case connectivity::boundary:
                zoneIsInternal[zonei] = false;
                break;
            case connectivity::internal:
                zoneIsInternal[zonei] = true;
                break;
            case connectivity::error:
                FatalErrorInFunction
                    << "Zone " << zoneNames[zonei]
                    << " contains both internal and boundary faces."
                    << " This is not allowed." << exit(FatalError);
                break;
        }
    }

    return zoneIsInternal;
}


// To combineReduce a labelList. Filters out duplicates.
struct uniqueEqOp
{
    void operator()(labelList& x, const labelList& y) const
    {
        if (x.empty())
        {
            if (y.size())
            {
                x = y;
            }
        }
        else
        {
            forAll(y, yi)
            {
                if (findIndex(x, y[yi]) == -1)
                {
                    label sz = x.size();
                    x.setSize(sz+1);
                    x[sz] = y[yi];
                }
            }
        }
    }
};


// Calculate global pp faces per pp edge.
labelListList globalEdgeFaces
(
    const polyMesh& mesh,
    const globalIndex& globalFaces,
    const primitiveFacePatch& pp,
    const labelList& ppMeshEdges
)
{
    // From mesh edge to global pp face labels.
    labelListList globalEdgeFaces(ppMeshEdges.size());

    const labelListList& edgeFaces = pp.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        // Store pp face and processor as unique tag.
        labelList& globalEFaces = globalEdgeFaces[edgeI];
        globalEFaces.setSize(eFaces.size());
        forAll(eFaces, i)
        {
            globalEFaces[i] = globalFaces.toGlobal(eFaces[i]);
        }
    }

    // Synchronise across coupled edges.
    syncTools::syncEdgeList
    (
        mesh,
        ppMeshEdges,
        globalEdgeFaces,
        uniqueEqOp(),
        labelList()             // null value
    );

    return globalEdgeFaces;
}


// Find a patch face that is not extruded. Return -1 if not found.
label findUncoveredPatchFace
(
    const fvMesh& mesh,
    const UIndirectList<label>& extrudeFaces, // mesh faces that are extruded
    const label meshEdgeI                     // mesh edge
)
{
    // Make set of extruded faces.
    labelHashSet extrudeFaceSet(extrudeFaces.size());
    forAll(extrudeFaces, i)
    {
        extrudeFaceSet.insert(extrudeFaces[i]);
    }

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];
    forAll(eFaces, i)
    {
        label facei = eFaces[i];
        label patchi = pbm.whichPatch(facei);

        if
        (
            patchi != -1
        && !pbm[patchi].coupled()
        && !extrudeFaceSet.found(facei)
        )
        {
            return facei;
        }
    }
    return -1;
}


// Same as findUncoveredPatchFace, except explicitly checks for cyclic faces
label findUncoveredCyclicPatchFace
(
    const fvMesh& mesh,
    const UIndirectList<label>& extrudeFaces, // mesh faces that are extruded
    const label meshEdgeI                     // mesh edge
)
{
    // Make set of extruded faces.
    labelHashSet extrudeFaceSet(extrudeFaces.size());
    forAll(extrudeFaces, i)
    {
        extrudeFaceSet.insert(extrudeFaces[i]);
    }

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];
    forAll(eFaces, i)
    {
        label facei = eFaces[i];
        label patchi = pbm.whichPatch(facei);

        if
        (
            patchi != -1
        &&  isA<cyclicPolyPatch>(pbm[patchi])
        && !extrudeFaceSet.found(facei)
        )
        {
            return facei;
        }
    }
    return -1;
}


// Add coupled patches into the mesh
void addCouplingPatches
(
    const fvMesh& mesh,
    const bool isShellMesh,
    const bool intrude,
    const word& regionName,
    const word& nbrRegionName,
    const wordList& zoneNames,
    const wordList& oppositeZoneNames,
    const boolList& zoneIsInternal,
    const dictionary& dict,
    DynamicList<polyPatch*>& newPatches,
    labelList& zoneTopPatch,
    labelList& zoneBottomPatch
)
{
    Pout<< "Adding coupling patches:" << nl << nl
        << "patchID\tpatch\ttype" << nl
        << "-------\t-----\t----"
        << endl;

    wordList patchNames
    (
        dict.lookupOrDefault("patchNames", wordList())
    );
    wordList regionPatchNames
    (
        dict.lookupOrDefault("regionPatchNames", wordList())
    );
    if (isShellMesh)
    {
        patchNames.swap(regionPatchNames);
    }

    const wordList patchTypes
    (
        isShellMesh
      ? dict.lookupOrDefault("regionPatchTypes", wordList())
      : dict.lookupOrDefault("patchTypes", wordList())
    );

    wordList oppositePatchNames
    (
        dict.lookupOrDefault("oppositePatchNames", wordList())
    );
    wordList regionOppositePatchNames
    (
        dict.lookupOrDefault("regionOppositePatchNames", wordList())
    );
    if (isShellMesh)
    {
        oppositePatchNames.swap(regionOppositePatchNames);
    }

    const wordList oppositePatchTypes
    (
        isShellMesh
      ? dict.lookupOrDefault("regionOppositePatchTypes", wordList())
      : dict.lookupOrDefault("oppositePatchType", wordList())
    );

    zoneTopPatch.setSize(zoneNames.size(), -1);
    zoneBottomPatch.setSize(zoneNames.size(), -1);

    dictionary patchDict;
    patchDict.add("neighbourRegion", nbrRegionName);

    label nOldPatches = newPatches.size();
    forAll(zoneNames, zonei)
    {
        const word patchNamePrefix =
            regionName + "_to_" + nbrRegionName + '_';

        const word nbrPatchNamePrefix =
            nbrRegionName + "_to_" + regionName + '_';

        word bottomPatchName;
        word bottomNbrPatchName;
        word topPatchName;
        word topNbrPatchName;

        if (zoneIsInternal[zonei])
        {
            bottomPatchName =
                patchNamePrefix + zoneNames[zonei] + "_bottom";
            bottomNbrPatchName =
                nbrPatchNamePrefix + zoneNames[zonei] + "_bottom";
            topPatchName =
                patchNamePrefix + zoneNames[zonei] + "_top";
            topNbrPatchName =
                nbrPatchNamePrefix + zoneNames[zonei] + "_top";
        }
        else if (!oppositeZoneNames[zonei].empty())
        {
            bottomPatchName = patchNamePrefix + zoneNames[zonei];
            bottomNbrPatchName = nbrPatchNamePrefix + zoneNames[zonei];
            topPatchName = patchNamePrefix + oppositeZoneNames[zonei];
            topNbrPatchName = nbrPatchNamePrefix + oppositeZoneNames[zonei];
        }
        else
        {
            bottomPatchName = patchNamePrefix + zoneNames[zonei];
            bottomNbrPatchName = nbrPatchNamePrefix + zoneNames[zonei];
            topPatchName = zoneNames[zonei] + "_top";
            topNbrPatchName = word::null;
        }

        if (patchNames.size())
        {
            bottomPatchName = patchNames[zonei];
        }
        if (regionPatchNames.size())
        {
            bottomNbrPatchName = regionPatchNames[zonei];
        }
        if (oppositePatchNames.size())
        {
            topPatchName = oppositePatchNames[zonei];
        }
        if (regionOppositePatchNames.size())
        {
            topNbrPatchName = regionOppositePatchNames[zonei];
        }

        const bool bothMapped =
            zoneIsInternal[zonei] || !oppositeZoneNames[zonei].empty();

        const bool bottomMapped = bothMapped || !(intrude && isShellMesh);
        const bool topMapped = bothMapped || intrude;

        const bool bottomExtruded = !isShellMesh && intrude;
        const bool topExtruded = !bottomExtruded;

        dictionary bottomPatchDict;
        if (bottomMapped)
        {
            bottomPatchDict = patchDict;
            bottomPatchDict.add
            (
                "neighbourPatch",
                !intrude ? bottomNbrPatchName : topNbrPatchName
            );
            if (bottomExtruded)
            {
                bottomPatchDict.add("isExtrudedRegion", isShellMesh);
            }
        }

        zoneBottomPatch[zonei] = addPatch
        (
            mesh.boundaryMesh(),
            bottomPatchName,
            patchTypes.size() ? patchTypes[zonei]
          : !bottomMapped ? polyPatch::typeName
          : !bottomExtruded
          ? mappedWallPolyPatch::typeName
          : mappedExtrudedWallPolyPatch::typeName,
            bottomPatchDict,
            newPatches
        );

        Pout<< zoneBottomPatch[zonei]
            << '\t' << newPatches[zoneBottomPatch[zonei]]->name()
            << '\t' << newPatches[zoneBottomPatch[zonei]]->type()
            << nl;

        if (!isShellMesh && !bothMapped) continue;

        dictionary topPatchDict;
        if (topMapped)
        {
            topPatchDict = patchDict;
            topPatchDict.add
            (
                "neighbourPatch",
                !intrude ? topNbrPatchName : bottomNbrPatchName
            );
            if (topExtruded)
            {
                topPatchDict.add("isExtrudedRegion", isShellMesh);
            }
        }

        zoneTopPatch[zonei] = addPatch
        (
            mesh.boundaryMesh(),
            topPatchName,
            oppositePatchTypes.size() ? oppositePatchTypes[zonei]
          : !topMapped ? polyPatch::typeName
          : !topExtruded
          ? mappedWallPolyPatch::typeName
          : mappedExtrudedWallPolyPatch::typeName,
            topPatchDict,
            newPatches
        );

        Pout<< zoneTopPatch[zonei]
            << '\t' << newPatches[zoneTopPatch[zonei]]->name()
            << '\t' << newPatches[zoneTopPatch[zonei]]->type()
            << nl;
    }

    Pout<< "Added " << newPatches.size()-nOldPatches
        << " inter-region patches." << nl
        << endl;
}


// Count the number of faces in patches that need to be created
labelList countExtrudePatches
(
    const fvMesh& mesh,
    const label nZones,
    const primitiveFacePatch& extrudePatch,
    const labelList& extrudeFaces,
    const labelList& extrudeFaceZoneIDs,
    const labelList& extrudeEdgeMeshEdges,
    const labelListList& extrudeEdgeGlobalFaces
)
{
    labelList zoneSideNFaces(nZones, 0);

    forAll(extrudePatch.edgeFaces(), edgeI)
    {
        const labelList& eFaces = extrudePatch.edgeFaces()[edgeI];

        if (eFaces.size() == 2)
        {
            // Internal edge
        }
        else if
        (
            eFaces.size() == 1
         && extrudeEdgeGlobalFaces[edgeI].size() == 2
        )
        {
            // Coupled edge
        }
        else
        {
            // Perimeter edge. Check whether we are on a mesh edge with
            // external patches. If so choose any uncovered one. If none found
            // put face in undetermined zone 'side' patch.
            const label facei = findUncoveredPatchFace
            (
                mesh,
                UIndirectList<label>(extrudeFaces, eFaces),
                extrudeEdgeMeshEdges[edgeI]
            );

            if (facei == -1)
            {
                forAll(extrudePatch.edgeFaces()[edgeI], i)
                {
                    const label facei = extrudePatch.edgeFaces()[edgeI][i];
                    zoneSideNFaces[extrudeFaceZoneIDs[facei]] ++;
                }
            }
        }
    }

    Pstream::listCombineGather(zoneSideNFaces, plusEqOp<label>());
    Pstream::listCombineScatter(zoneSideNFaces);

    return zoneSideNFaces;
}


// Add side patches into the mesh
labelList addZoneSidePatches
(
    const fvMesh& mesh,
    const wordList& zoneNames,
    const labelList& zoneSideNFaces,
    const word& oneDPolyPatchType,
    DynamicList<polyPatch*>& newPatches
)
{
    Pout<< "Adding patches for sides on zones:" << nl << nl
        << "patchID\tpatch" << nl << "-------\t-----" << endl;

    labelList zoneSidePatches(zoneNames.size(), -labelMax);

    const label nOldPatches = newPatches.size();

    if (oneDPolyPatchType != word::null)
    {
        forAll(zoneSideNFaces, zoneI)
        {
            word patchName;

            if (oneDPolyPatchType == "empty")
            {
                patchName = "oneDEmptyPatch";
                zoneSidePatches[zoneI] = addPatch
                (
                    mesh.boundaryMesh(),
                    patchName,
                    emptyPolyPatch::typeName,
                    dictionary(),
                    newPatches
                );
            }
            else if (oneDPolyPatchType == "wedge")
            {
                patchName = "oneDWedgePatch";
                zoneSidePatches[zoneI] = addPatch
                (
                    mesh.boundaryMesh(),
                    patchName,
                    wedgePolyPatch::typeName,
                    dictionary(),
                    newPatches
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Type " << oneDPolyPatchType << " does not exist "
                    << exit(FatalError);
            }

            Pout<< zoneSidePatches[zoneI] << '\t' << patchName << nl;
        }
    }
    else
    {
        forAll(zoneSideNFaces, zoneI)
        {
            if (zoneSideNFaces[zoneI] > 0)
            {
                word patchName = zoneNames[zoneI] + "_" + "side";

                zoneSidePatches[zoneI] = addPatch
                (
                    mesh.boundaryMesh(),
                    patchName,
                    polyPatch::typeName,
                    dictionary(),
                    newPatches
                );

                Pout<< zoneSidePatches[zoneI] << '\t' << patchName << nl;
            }
        }
    }

    Pout<< "Added " << newPatches.size() - nOldPatches << " zone-side patches."
        << nl << endl;

    return zoneSidePatches;
}


// Add any additional coupled side patches that might be necessary. Return
// correspondence between extruded edges and their side patches.
labelList addExtrudeEdgeSidePatches
(
    const fvMesh& mesh,
    const primitiveFacePatch& extrudePatch,
    const labelList& extrudeFaces,
    const labelList& extrudeEdgeMeshEdges,
    const distributionMap& extrudeEdgeFacesMap,
    const labelListList& extrudeEdgeGlobalFaces,
    DynamicList<polyPatch*>& newPatches
)
{
    // Get procID in extrudeEdgeGlobalFaces order
    labelList procID(extrudeEdgeGlobalFaces.size(), Pstream::myProcNo());
    extrudeEdgeFacesMap.distribute(procID);

    // Calculate opposite processor for coupled edges (only if shared by
    // two procs). Note: Could have saved original globalEdgeFaces structure.
    labelList minProcID(extrudeEdgeGlobalFaces.size(), labelMax);
    labelList maxProcID(extrudeEdgeGlobalFaces.size(), labelMin);
    forAll(extrudeEdgeGlobalFaces, edgeI)
    {
        const labelList& eFaces = extrudeEdgeGlobalFaces[edgeI];
        if (eFaces.size())
        {
            forAll(eFaces, i)
            {
                label proci = procID[eFaces[i]];
                minProcID[edgeI] = min(minProcID[edgeI], proci);
                maxProcID[edgeI] = max(maxProcID[edgeI], proci);
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh,
        extrudeEdgeMeshEdges,
        minProcID,
        minEqOp<label>(),
        labelMax        // null value
    );
    syncTools::syncEdgeList
    (
        mesh,
        extrudeEdgeMeshEdges,
        maxProcID,
        maxEqOp<label>(),
        labelMin        // null value
    );

    Pout<< "Adding processor or cyclic patches:" << nl << nl
        << "patchID\tpatch" << nl << "-------\t-----" << endl;

    const label nOldPatches = newPatches.size();

    labelList extrudeEdgeSidePatches(extrudePatch.edgeFaces().size(), -1);
    forAll(extrudePatch.edgeFaces(), edgeI)
    {
        const labelList& eFaces = extrudePatch.edgeFaces()[edgeI];

        if
        (
            eFaces.size() == 1
         && extrudeEdgeGlobalFaces[edgeI].size() == 2
        )
        {
            // Coupled boundary edge. Find matching patch.
            label nbrProci = minProcID[edgeI];
            if (nbrProci == Pstream::myProcNo())
            {
                nbrProci = maxProcID[edgeI];
            }

            if (nbrProci == Pstream::myProcNo())
            {
                // Cyclic patch since both procs the same. This cyclic should
                // already exist in newPatches so no adding necessary.

                const label facei = findUncoveredCyclicPatchFace
                (
                    mesh,
                    UIndirectList<label>(extrudeFaces, eFaces),
                    extrudeEdgeMeshEdges[edgeI]
                );

                if (facei != -1)
                {
                    const polyBoundaryMesh& patches = mesh.boundaryMesh();

                    const label newPatchi = findPatchID
                    (
                        newPatches,
                        patches[patches.whichPatch(facei)].name()
                    );

                    extrudeEdgeSidePatches[edgeI] = newPatchi;
                }
                else
                {
                    FatalErrorInFunction
                        << "Unable to determine coupled patch addressing"
                        << abort(FatalError);
                }
            }
            else
            {
                // Processor patch
                const word name =
                    processorPolyPatch::newName(Pstream::myProcNo(), nbrProci);

                extrudeEdgeSidePatches[edgeI] = findPatchID(newPatches, name);

                if (extrudeEdgeSidePatches[edgeI] == -1)
                {
                    dictionary patchDict;
                    patchDict.add("myProcNo", Pstream::myProcNo());
                    patchDict.add("neighbProcNo", nbrProci);

                    extrudeEdgeSidePatches[edgeI] = addPatch
                    (
                        mesh.boundaryMesh(),
                        name,
                        processorPolyPatch::typeName,
                        patchDict,
                        newPatches
                    );

                    Pout<< extrudeEdgeSidePatches[edgeI] << '\t' << name << nl;
                }
            }
        }
    }

    Pout<< "Added " << newPatches.size() - nOldPatches
        << " coupled patches." << nl << endl;

    return extrudeEdgeSidePatches;
}


int main(int argc, char *argv[])
{
    argList::addNote("Create region mesh by extruding a faceZone or faceSet");

    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    if (mesh.boundaryMesh().checkParallelSync(true))
    {
        List<wordList> allNames(Pstream::nProcs());
        allNames[Pstream::myProcNo()] = mesh.boundaryMesh().names();
        Pstream::gatherList(allNames);
        Pstream::scatterList(allNames);

        FatalErrorInFunction
            << "Patches are not synchronised on all processors."
            << " Per processor patches " << allNames
            << exit(FatalError);
    }

    bool overwrite = args.optionFound("overwrite");

    const dictionary dict(systemDict("extrudeToRegionMeshDict", args, mesh));

    // Region to create by extrusion
    const word shellRegionName(dict.lookup("region"));

    // Should the extruded region overlap the existing region, i.e. "intrude"?
    const Switch intrude(dict.lookupOrDefault("intrude", false));

    if (shellRegionName == regionName)
    {
        FatalIOErrorIn(args.executable().c_str(), dict)
            << "Cannot extrude into same region as mesh." << endl
            << "Mesh region : " << regionName << endl
            << "Shell region : " << shellRegionName
            << exit(FatalIOError);
    }

    // Select faces to extrude
    enum class zoneSourceType
    {
        zone,
        set,
        patch
    };

    static const wordList zoneSourceTypeNames
    {
        "faceZone",
        "faceSet",
        "patch"
    };

    static const wordList zoneSourcesTypeNames
    {
        "faceZones",
        "faceSets",
        "patches"
    };

    wordList zoneNames;
    wordList oppositeZoneNames;
    List<zoneSourceType> zoneSourceTypes;

    auto lookupZones = [&](const zoneSourceType& type)
    {
        const word& keyword = zoneSourcesTypeNames[unsigned(type)];

        if (dict.found(keyword))
        {
            zoneNames.append(dict.lookup<wordList>(keyword));

            oppositeZoneNames.append
            (
                dict.lookupOrDefaultBackwardsCompatible
                (
                    {
                        "opposite" + keyword.capitalise(),
                        keyword + "Shadow"
                    },
                    wordList
                    (
                        zoneNames.size() - oppositeZoneNames.size(),
                        word::null
                    )
                )
            );

            zoneSourceTypes.setSize(zoneNames.size(), type);
        }
    };
    lookupZones(zoneSourceType::zone);
    lookupZones(zoneSourceType::set);
    lookupZones(zoneSourceType::patch);

    Info<< nl << "Extruding:" << nl << incrIndent;
    forAll(zoneNames, zonei)
    {
        const unsigned typei = unsigned(zoneSourceTypes[zonei]);

        if (oppositeZoneNames[zonei].empty())
        {
            Info<< indent << "From " << zoneSourceTypeNames[typei] << " \""
                << zoneNames[zonei] << "\"" << nl;
        }
        else
        {
            Info<< indent << "Between " << zoneSourcesTypeNames[typei] << " \""
                << zoneNames[zonei] << "\" and \"" << oppositeZoneNames[zonei]
                << "\"" << nl;
        }
    }
    Info<< endl << decrIndent;

    // One-dimensional extrusion settings
    const Switch oneD(dict.lookupOrDefault("oneD", false));
    Switch oneDNonManifoldEdges(false);
    word oneDPatchType(emptyPolyPatch::typeName);
    if (oneD)
    {
        oneDNonManifoldEdges = dict.lookupOrDefault("nonManifold", false);
        dict.lookup("oneDPolyPatchType") >> oneDPatchType;
    }

    if (oneD)
    {
        if (oneDNonManifoldEdges)
        {
            Info<< "Extruding as 1D columns with sides in patch type "
                << oneDPatchType
                << " and connected points (except on non-manifold areas)."
                << endl;
        }
        else
        {
            Info<< "Extruding as 1D columns with sides in patch type "
                << oneDPatchType
                << " and duplicated points (overlapping volumes)."
                << endl;
        }
    }

    // Construct the point generator
    autoPtr<extrudeModel> model(extrudeModel::New(dict));

    // Change the primary mesh?
    const Switch adaptMesh(dict.lookup("adaptMesh"));


    // Determine output instance
    word meshInstance;
    if (!overwrite)
    {
        runTime++;
        meshInstance = runTime.name();
    }
    else
    {
        meshInstance = mesh.pointsInstance();
    }
    Info<< "Writing meshes to " << meshInstance << nl << endl;


    // Map from extrude zone to mesh zone, or -1 if not a mesh zone
    labelList zoneMeshZoneID(zoneNames.size(), -1);
    labelList oppositeZoneMeshZoneID(zoneNames.size(), -1);
    forAll(zoneNames, zonei)
    {
        if (zoneSourceTypes[zonei] != zoneSourceType::zone) continue;

        zoneMeshZoneID[zonei] =
            mesh.faceZones().findZoneID(zoneNames[zonei]);

        if (zoneMeshZoneID[zonei] == -1)
        {
            FatalIOErrorIn(args.executable().c_str(), dict)
                << "Cannot find zone " << zoneNames[zonei]
                << endl << "Valid zones are " << mesh.faceZones().names()
                << exit(FatalIOError);
        }

        if (!oppositeZoneNames[zonei].empty())
        {
            oppositeZoneMeshZoneID[zonei] =
                mesh.faceZones().findZoneID(oppositeZoneNames[zonei]);

            if (oppositeZoneMeshZoneID[zonei] == -1)
            {
                FatalIOErrorIn(args.executable().c_str(), dict)
                    << "Cannot find opposite zone " << oppositeZoneNames[zonei]
                    << endl << "Valid zones are " << mesh.faceZones().names()
                    << exit(FatalIOError);
            }
        }
    }


    // Extract faces to extrude
    labelList extrudeFaces, oppositeExtrudeFaces;
    labelList extrudeFaceZoneIDs, oppositeExtrudeFaceZoneIDs;
    boolList extrudeFaceFlips, oppositeExtrudeFaceFlips;
    {
        // Load any faceSets that we need
        PtrList<faceSet> zoneSets(zoneNames.size());
        PtrList<faceSet> oppositeZoneSets(zoneNames.size());
        forAll(zoneNames, zonei)
        {
            if (zoneSourceTypes[zonei] == zoneSourceType::set)
            {
                zoneSets.set(zonei, new faceSet(mesh, zoneNames[zonei]));
                if (!oppositeZoneNames.empty())
                {
                    oppositeZoneSets.set
                    (
                       zonei,
                       new faceSet(mesh, oppositeZoneNames[zonei])
                    );
                }
            }
        }

        // Create dynamic face lists
        DynamicList<label> facesDyn, oppositeFacesDyn;
        DynamicList<label> zoneIDsDyn, oppositeZoneIDsDyn;
        DynamicList<bool> flipsDyn, oppositeFlipsDyn;
        forAll(zoneNames, zonei)
        {
            switch (zoneSourceTypes[zonei])
            {
                case zoneSourceType::zone:
                {
                    const faceZone& fz =
                        mesh.faceZones()[zoneMeshZoneID[zonei]];
                    facesDyn.append(fz);
                    zoneIDsDyn.append(labelList(fz.size(), zonei));
                    flipsDyn.append(fz.flipMap());

                    if (!oppositeZoneNames[zonei].empty())
                    {
                        const faceZone& sfz =
                            mesh.faceZones()[oppositeZoneMeshZoneID[zonei]];
                        if (sfz.size() != fz.size())
                        {
                            FatalIOErrorIn(args.executable().c_str(), dict)
                                << "Opposite zone " << oppositeZoneNames[zonei]
                                << "is a different size from it's "
                                << "corresponding zone " << zoneNames[zonei]
                                << exit(FatalIOError);
                        }
                        oppositeFacesDyn.append(sfz);
                        oppositeZoneIDsDyn.append(labelList(sfz.size(), zonei));
                        oppositeFlipsDyn.append(sfz.flipMap());
                    }
                    else
                    {
                        oppositeFacesDyn.append(labelList(fz.size(), -1));
                        oppositeZoneIDsDyn.append(labelList(fz.size(), -1));
                        oppositeFlipsDyn.append(boolList(fz.size(), false));
                    }
                    break;
                }
                case zoneSourceType::set:
                {
                    const faceSet& fs = zoneSets[zonei];
                    facesDyn.append(fs.toc());
                    zoneIDsDyn.append(labelList(fs.size(), zonei));
                    flipsDyn.append(boolList(fs.size(), false));

                    if (!oppositeZoneNames[zonei].empty())
                    {
                        const faceSet& sfs = oppositeZoneSets[zonei];
                        if (sfs.size() != fs.size())
                        {
                            FatalIOErrorIn(args.executable().c_str(), dict)
                                << "Opposite set " << oppositeZoneNames[zonei]
                                << "is a different size from it's "
                                << "corresponding zone " << zoneNames[zonei]
                                << exit(FatalIOError);
                        }
                        oppositeFacesDyn.append(sfs.toc());
                        oppositeZoneIDsDyn.append(labelList(sfs.size(), zonei));
                        oppositeFlipsDyn.append(boolList(sfs.size(), false));
                    }
                    else
                    {
                        oppositeFacesDyn.append(labelList(fs.size(), -1));
                        oppositeZoneIDsDyn.append(labelList(fs.size(), -1));
                        oppositeFlipsDyn.append(boolList(fs.size(), false));
                    }
                    break;
                }
                case zoneSourceType::patch:
                {
                    const polyPatch& pp =
                        mesh.boundaryMesh()[zoneNames[zonei]];
                    facesDyn.append(pp.start() + identityMap(pp.size()));
                    zoneIDsDyn.append(labelList(pp.size(), zonei));
                    flipsDyn.append(boolList(pp.size(), false));

                    if (!oppositeZoneNames[zonei].empty())
                    {
                        const polyPatch& spp =
                            mesh.boundaryMesh()[oppositeZoneNames[zonei]];
                        if (spp.size() != pp.size())
                        {
                            FatalIOErrorIn(args.executable().c_str(), dict)
                                << "Opposite patch " << oppositeZoneNames[zonei]
                                << "is a different size from it's "
                                << "corresponding zone " << zoneNames[zonei]
                                << exit(FatalIOError);
                        }
                        oppositeFacesDyn.append
                        (
                            spp.start() + identityMap(spp.size())
                        );
                        oppositeZoneIDsDyn.append(labelList(spp.size(), zonei));
                        oppositeFlipsDyn.append(boolList(spp.size(), false));
                    }
                    else
                    {
                        oppositeFacesDyn.append(labelList(pp.size(), -1));
                        oppositeZoneIDsDyn.append(labelList(pp.size(), -1));
                        oppositeFlipsDyn.append(boolList(pp.size(), false));
                    }
                    break;
                }
            }
        }

        // Transfer to non-dynamic storage
        extrudeFaces.transfer(facesDyn);
        extrudeFaceZoneIDs.transfer(zoneIDsDyn);
        extrudeFaceFlips.transfer(flipsDyn);
        oppositeExtrudeFaces.transfer(oppositeFacesDyn);
        oppositeExtrudeFaceZoneIDs.transfer(oppositeZoneIDsDyn);
        oppositeExtrudeFaceFlips.transfer(oppositeFlipsDyn);
    }

    faceList extrudeFaceList(UIndirectList<face>(mesh.faces(), extrudeFaces));

    forAll(extrudeFaceList, facei)
    {
        if (intrude != extrudeFaceFlips[facei])
        {
            extrudeFaceList[facei].flip();
        }
    }

    // Create a primitive patch of the extruded faces
    const primitiveFacePatch extrudePatch
    (
        extrudeFaceList,
        mesh.points()
    );

    // Check zone either all internal or all external faces
    const boolList zoneIsInternal
    (
        calcZoneIsInternal(mesh, zoneNames, extrudeFaces, extrudeFaceZoneIDs)
    );


    Info<< "extrudePatch :"
        << " faces:" << returnReduce(extrudePatch.size(), sumOp<label>())
        << " points:" << returnReduce(extrudePatch.nPoints(), sumOp<label>())
        << " edges:" << returnReduce(extrudePatch.nEdges(), sumOp<label>())
        << nl << endl;


    // Determine per-extrude-edge info
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Corresponding mesh edges
    const labelList extrudeEdgeMeshEdges
    (
        extrudePatch.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    const globalIndex globalExtrudeFaces(extrudePatch.size());

    // Global pp faces per pp edge.
    labelListList extrudeEdgeGlobalFaces
    (
        globalEdgeFaces
        (
            mesh,
            globalExtrudeFaces,
            extrudePatch,
            extrudeEdgeMeshEdges
        )
    );
    List<Map<label>> compactMap;
    const distributionMap extrudeEdgeFacesMap
    (
        globalExtrudeFaces,
        extrudeEdgeGlobalFaces,
        compactMap
    );


    // Copy all non-local patches since these are used on boundary edges of
    // the extrusion
    DynamicList<polyPatch*> regionPatches(mesh.boundaryMesh().size());
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (!isA<processorPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            regionPatches.append
            (
                mesh.boundaryMesh()[patchi].clone
                (
                    mesh.boundaryMesh(),
                    regionPatches.size(),
                    0,              // size
                    0               // start
                ).ptr()
            );
        }
    }


    // Add interface patches
    // ~~~~~~~~~~~~~~~~~~~~~

    // From zone to interface patch (region side)
    labelList zoneTopPatch, zoneBottomPatch;
    addCouplingPatches
    (
        mesh,
        true,
        intrude,
        shellRegionName,
        regionName,
        zoneNames,
        oppositeZoneNames,
        zoneIsInternal,
        dict,
        regionPatches,
        zoneTopPatch,
        zoneBottomPatch
    );

    // From zone to interface patch (mesh side)
    labelList interMeshTopPatch;
    labelList interMeshBottomPatch;
    if (adaptMesh)
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Clone existing non-processor patches
        DynamicList<polyPatch*> newPatches(patches.size());
        forAll(patches, patchi)
        {
            if (!isA<processorPolyPatch>(patches[patchi]))
            {
                newPatches.append(patches[patchi].clone(patches).ptr());
            }
        }

        // Add new patches
        addCouplingPatches
        (
            mesh,
            false,
            intrude,
            regionName,
            shellRegionName,
            zoneNames,
            oppositeZoneNames,
            zoneIsInternal,
            dict,
            newPatches,
            interMeshTopPatch,
            interMeshBottomPatch
        );

        // Clone existing processor patches
        forAll(patches, patchi)
        {
            if (isA<processorPolyPatch>(patches[patchi]))
            {
                newPatches.append
                (
                    patches[patchi].clone
                    (
                        patches,
                        newPatches.size(),
                        patches[patchi].size(),
                        patches[patchi].start()
                    ).ptr()
                );
            }
        }

        // Add to mesh
        mesh.clearOut();
        mesh.removeFvBoundary();
        mesh.addFvPatches(newPatches, true);

        // Note: from this point on mesh patches differs from regionPatches
    }


    // Patch per extruded face
    labelList extrudeFaceTopPatchID(extrudePatch.size());
    labelList extrudeFaceBottomPatchID(extrudePatch.size());
    forAll(extrudeFaceZoneIDs, facei)
    {
        extrudeFaceTopPatchID[facei] =
            zoneTopPatch[extrudeFaceZoneIDs[facei]];
        extrudeFaceBottomPatchID[facei] =
            zoneBottomPatch[extrudeFaceZoneIDs[facei]];
    }


    // Count how many patches on special edges of extrudePatch are necessary
    const labelList zoneSideNFaces
    (
        countExtrudePatches
        (
            mesh,
            zoneNames.size(),
            extrudePatch,           // patch
            extrudeFaces,           // mesh face per patch face
            extrudeFaceZoneIDs,     // ...
            extrudeEdgeMeshEdges,   // mesh edge per patch edge
            extrudeEdgeGlobalFaces  // global indexing per patch edge
        )
    );


    // Add the zone-side patches.
    const labelList zoneSidePatches
    (
        addZoneSidePatches
        (
            mesh,
            zoneNames,
            zoneSideNFaces,
            (oneD ? oneDPatchType : word::null),
            regionPatches
        )
    );


    // Sets extrudeEdgeSidePatches[edgei] to interprocessor patch. Adds any
    // interprocessor or cyclic patches if necessary.
    const labelList extrudeEdgeSidePatches
    (
        addExtrudeEdgeSidePatches
        (
            mesh,
            extrudePatch,
            extrudeFaces,
            extrudeEdgeMeshEdges,
            extrudeEdgeFacesMap,
            extrudeEdgeGlobalFaces,
            regionPatches
        )
    );


    // Set patches to use for edges to be extruded into boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // In order of edgeFaces: per edge, per originating face the
    // patch to use for the side face (from the extruded edge).
    // If empty size create an internal face.
    labelListList extrudeEdgePatches(extrudePatch.nEdges());

    // Is edge a non-manifold edge
    PackedBoolList nonManifoldEdge(extrudePatch.nEdges(), false);

    // Note: logic has to be same as in countExtrudePatches.
    forAll(extrudePatch.edgeFaces(), edgeI)
    {
        const labelList& eFaces = extrudePatch.edgeFaces()[edgeI];
        labelList& ePatches = extrudeEdgePatches[edgeI];

        if (oneD)
        {
            ePatches.setSize(eFaces.size());
            forAll(eFaces, i)
            {
                ePatches[i] =
                    zoneSidePatches[extrudeFaceZoneIDs[eFaces[i]]];
            }

            if (oneDNonManifoldEdges)
            {
                if (eFaces.size() != 2)
                {
                    nonManifoldEdge[edgeI] = true;
                }
            }
            else
            {
                nonManifoldEdge[edgeI] = true;
            }
        }
        else if (eFaces.size() == 2)
        {
            // Internal edge
        }
        else if (extrudeEdgeSidePatches[edgeI] != -1)
        {
            // Coupled edge
            ePatches.setSize(eFaces.size());
            forAll(eFaces, i)
            {
                ePatches[i] = extrudeEdgeSidePatches[edgeI];
            }
        }
        else
        {
            // Perimeter edge
            const label facei = findUncoveredPatchFace
            (
                mesh,
                UIndirectList<label>(extrudeFaces, eFaces),
                extrudeEdgeMeshEdges[edgeI]
            );

            if (facei != -1)
            {
                const label newPatchi =
                    findPatchID
                    (
                        regionPatches,
                        mesh.boundaryMesh()
                        [
                            mesh.boundaryMesh().whichPatch(facei)
                        ].name()
                    );
                ePatches.setSize(eFaces.size(), newPatchi);
            }
            else
            {
                ePatches.setSize(eFaces.size());
                forAll(eFaces, i)
                {
                    ePatches[i] =
                        zoneSidePatches[extrudeFaceZoneIDs[eFaces[i]]];
                }
            }

            nonManifoldEdge[edgeI] = true;
        }
    }



    // Assign point regions
    // ~~~~~~~~~~~~~~~~~~~~

    // Per face, per point the region number.
    faceList pointGlobalRegions;
    faceList pointLocalRegions;
    labelList localToGlobalRegion;

    createShellMesh::calcPointRegions
    (
        mesh.globalData(),
        extrudePatch,
        nonManifoldEdge,
        false,              // keep cyclic separated regions apart
        pointGlobalRegions,
        pointLocalRegions,
        localToGlobalRegion
    );

    // Per local region an originating point
    labelList localRegionPoints(localToGlobalRegion.size());
    forAll(pointLocalRegions, facei)
    {
        const face& f = extrudePatch.localFaces()[facei];
        const face& pRegions = pointLocalRegions[facei];
        forAll(pRegions, fp)
        {
            localRegionPoints[pRegions[fp]] = f[fp];
        }
    }

    // Calculate region normals by reducing local region normals
    pointField localRegionNormals(localToGlobalRegion.size());
    {
        pointField localSum(localToGlobalRegion.size(), Zero);

        forAll(pointLocalRegions, facei)
        {
            const face& pRegions = pointLocalRegions[facei];
            forAll(pRegions, fp)
            {
                const label localRegionI = pRegions[fp];

                // Add a small amount of the face-centre-to-point vector in
                // order to stabilise the computation of normals on the edges
                // of baffles
                localSum[localRegionI] +=
                    rootSmall
                   *(
                       extrudePatch.points()[extrudePatch[facei][fp]]
                     - extrudePatch.faceCentres()[facei]
                    )
                  + (1 - rootSmall)
                   *(
                        extrudePatch.faceNormals()[facei]
                    );
            }
        }

        Map<point> globalSum(2*localToGlobalRegion.size());

        forAll(localSum, localRegionI)
        {
            label globalRegionI = localToGlobalRegion[localRegionI];
            globalSum.insert(globalRegionI, localSum[localRegionI]);
        }

        // Reduce
        Pstream::mapCombineGather(globalSum, plusEqOp<point>());
        Pstream::mapCombineScatter(globalSum);

        forAll(localToGlobalRegion, localRegionI)
        {
            label globalRegionI = localToGlobalRegion[localRegionI];
            localRegionNormals[localRegionI] = globalSum[globalRegionI];
        }
        localRegionNormals /= mag(localRegionNormals) ;
    }


    // Use model to create displacements of first layer
    vectorField firstDisp(localRegionNormals.size());
    forAll(firstDisp, regionI)
    {
        const point& regionPt = extrudePatch.points()
        [
            extrudePatch.meshPoints()
            [
                localRegionPoints[regionI]
            ]
        ];
        const vector& n = localRegionNormals[regionI];
        firstDisp[regionI] = model()(regionPt, n, 1) - regionPt;
    }


    // Create a new mesh
    // ~~~~~~~~~~~~~~~~~

    fvMesh regionMesh
    (
        IOobject
        (
            shellRegionName,
            meshInstance,
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        pointField(),
        faceList(),
        labelList(),
        labelList(),
        false
    );

    // Add the new patches
    forAll(regionPatches, patchi)
    {
        polyPatch* ppPtr = regionPatches[patchi];
        regionPatches[patchi] = ppPtr->clone(regionMesh.boundaryMesh()).ptr();
        delete ppPtr;
    }
    regionMesh.clearOut();
    regionMesh.removeFvBoundary();
    regionMesh.addFvPatches(regionPatches, true);

    // Extrude
    {
        polyTopoChange meshMod(regionPatches.size());

        createShellMesh extruder
        (
            extrudePatch,
            pointLocalRegions,
            localRegionPoints
        );

        extruder.setRefinement
        (
            firstDisp,                              // first displacement
            model().expansionRatio(),
            model().nLayers(),                      // nLayers
            extrudeFaceTopPatchID,
            extrudeFaceBottomPatchID,
            extrudeEdgePatches,
            meshMod
        );

        // Enforce actual point positions according to extrudeModel (model), as
        // the extruder only does a fixed expansionRatio. The regionPoints and
        // nLayers are looped in the same way as in createShellMesh.
        DynamicList<point>& newPoints =
            const_cast<DynamicList<point>&>(meshMod.points());
        label meshPointi = extrudePatch.localPoints().size();
        forAll(localRegionPoints, regioni)
        {
            const label pointi = localRegionPoints[regioni];
            const point& pt = extrudePatch.localPoints()[pointi];
            const vector& n = localRegionNormals[regioni];

            for (label layeri = 1; layeri <= model().nLayers(); layeri++)
            {
                newPoints[meshPointi++] = model()(pt, n, layeri);
            }
        }

        meshMod.changeMesh(regionMesh, false);
    }

    // Set region mesh instance and write options
    regionMesh.setInstance(meshInstance);

    // Remove any unused patches
    deleteEmptyPatches(regionMesh);


    Info<< "Writing mesh " << regionMesh.name() << " to "
        << regionMesh.facesInstance() << nl << endl;
    regionMesh.write();


    // Insert baffles into original mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<polyTopoChangeMap> addBafflesMap;

    if (adaptMesh)
    {
        polyTopoChange meshMod(mesh);

        // Modify faces to be in bottom (= always coupled) patch
        forAll(extrudeFaces, facei)
        {
            const label meshFacei = extrudeFaces[facei];
            const label zonei = extrudeFaceZoneIDs[facei];
            const bool flip = extrudeFaceFlips[facei];
            const face& f = mesh.faces()[meshFacei];

            if (!flip)
            {
                meshMod.modifyFace
                (
                    f,                          // modified face
                    meshFacei,                  // label of face being modified
                    mesh.faceOwner()[meshFacei],// owner
                    -1,                         // neighbour
                    false,                      // face flip
                    interMeshBottomPatch[zonei],// patch for face
                    zoneMeshZoneID[zonei],      // zone for face
                    false                       // face flip in zone
                );
            }
            else if (mesh.isInternalFace(meshFacei))
            {
                meshMod.modifyFace
                (
                    f.reverseFace(),                // modified face
                    meshFacei,                      // label of modified face
                    mesh.faceNeighbour()[meshFacei],// owner
                    -1,                             // neighbour
                    true,                           // face flip
                    interMeshBottomPatch[zonei],    // patch for face
                    zoneMeshZoneID[zonei],          // zone for face
                    true                            // face flip in zone
                );
            }
        }

        forAll(extrudeFaces, facei)
        {
            if (oppositeExtrudeFaces[facei] != -1)
            {
                const label meshFacei = oppositeExtrudeFaces[facei];
                const label zonei = oppositeExtrudeFaceZoneIDs[facei];
                const bool flip = oppositeExtrudeFaceFlips[facei];
                const face& f = mesh.faces()[meshFacei];

                if (!flip)
                {
                    meshMod.modifyFace
                    (
                        f,                          // modified face
                        meshFacei,                  // face being modified
                        mesh.faceOwner()[meshFacei],// owner
                        -1,                         // neighbour
                        false,                      // face flip
                        interMeshTopPatch[zonei],   // patch for face
                        zoneMeshZoneID[zonei],      // zone for face
                        false                       // face flip in zone
                    );
                }
                else if (mesh.isInternalFace(meshFacei))
                {
                    meshMod.modifyFace
                    (
                        f.reverseFace(),                // modified face
                        meshFacei,                      // label modified face
                        mesh.faceNeighbour()[meshFacei],// owner
                        -1,                             // neighbour
                        true,                           // face flip
                        interMeshTopPatch[zonei],       // patch for face
                        zoneMeshZoneID[zonei],          // zone for face
                        true                            // face flip in zone
                    );
                }
            }
            else
            {
                const label meshFacei = extrudeFaces[facei];
                const label zonei = extrudeFaceZoneIDs[facei];
                const bool flip = extrudeFaceFlips[facei];
                const face& f = mesh.faces()[meshFacei];

                if (!flip)
                {
                    if (mesh.isInternalFace(meshFacei))
                    {
                        meshMod.addFace
                        (
                            f.reverseFace(),                // modified face
                            mesh.faceNeighbour()[meshFacei],// owner
                            -1,                             // neighbour
                            -1,                             // master point
                            -1,                             // master edge
                            meshFacei,                      // master face
                            true,                           // flip flux
                            interMeshTopPatch[zonei],       // patch for face
                            -1,                             // zone for face
                            false                           // face flip in zone
                        );
                    }
                }
                else
                {
                    meshMod.addFace
                    (
                        f,                              // face
                        mesh.faceOwner()[meshFacei],    // owner
                        -1,                             // neighbour
                        -1,                             // master point
                        -1,                             // master edge
                        meshFacei,                      // master face
                        false,                          // flip flux
                        interMeshTopPatch[zonei],       // patch for face
                        -1,                             // zone for face
                        false                           // zone flip
                    );
                }
            }
        }

        // Change the mesh. Change points directly (no inflation).
        addBafflesMap = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.topoChange(addBafflesMap);

        // Move mesh (since morphing might not do this)
        if (addBafflesMap().hasMotionPoints())
        {
            mesh.setPoints(addBafflesMap().preMotionPoints());
        }

        mesh.setInstance(meshInstance);

        // Remove any unused patches
        deleteEmptyPatches(mesh);

        Info<< "Writing mesh " << mesh.name() << " to "
            << mesh.facesInstance() << nl << endl;
        mesh.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
