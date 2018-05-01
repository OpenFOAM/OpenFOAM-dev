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
    - extrude edges in between different faceZones as a
      (nonuniformTransform)cyclic \<zoneA\>_\<zoneB\>
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
#include "fvMesh.H"
#include "polyTopoChange.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mappedWallPolyPatch.H"
#include "createShellMesh.H"
#include "syncTools.H"
#include "cyclicPolyPatch.H"
#include "wedgePolyPatch.H"
#include "nonuniformTransformCyclicPolyPatch.H"
#include "extrudeModel.H"
#include "globalIndex.H"
#include "faceSet.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
//#include "ReadFields.H"
#include "fvMeshTools.H"
#include "OBJstream.H"
#include "PatchTools.H"

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


template<class PatchType>
label addPatch
(
    const polyBoundaryMesh& patches,
    const word& patchName,
    DynamicList<polyPatch*>& newPatches
)
{
    label patchi = findPatchID(newPatches, patchName);

    if (patchi != -1)
    {
        if (isA<PatchType>(*newPatches[patchi]))
        {
            // Already there
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


    newPatches.append
    (
        polyPatch::New
        (
            PatchType::typeName,
            patchName,
            0,                          // size
            startFacei,                 // nFaces
            patchi,
            patches
        ).ptr()
    );

    return patchi;
}


template<class PatchType>
label addPatch
(
    const polyBoundaryMesh& patches,
    const word& patchName,
    const dictionary& dict,
    DynamicList<polyPatch*>& newPatches
)
{
    label patchi = findPatchID(newPatches, patchName);

    if (patchi != -1)
    {
        if (isA<PatchType>(*newPatches[patchi]))
        {
            // Already there
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
    patchDict.set("type", PatchType::typeName);
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


void createDummyFvMeshFiles(const polyMesh& mesh, const word& regionName)
{
    // Create dummy system/fv*
    {
        IOobject io
        (
            "fvSchemes",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        Info<< "Testing:" << io.objectPath() << endl;

        if (!io.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(true))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
}


// Check zone either all internal or all external faces
void checkZoneInside
(
    const polyMesh& mesh,
    const wordList& zoneNames,
    const labelList& zoneID,
    const labelList& extrudeMeshFaces,
    const boolList& isInternal
)
{
    forAll(zoneNames, i)
    {
        if (isInternal[i])
        {
            Info<< "Zone " << zoneNames[i] << " has internal faces" << endl;
        }
        else
        {
            Info<< "Zone " << zoneNames[i] << " has boundary faces" << endl;
        }
    }

    forAll(extrudeMeshFaces, i)
    {
        label facei = extrudeMeshFaces[i];
        label zoneI = zoneID[i];
        if (isInternal[zoneI] != mesh.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Zone " << zoneNames[zoneI]
                << " is not consistently all internal or all boundary faces."
                << " Face " << facei << " at " << mesh.faceCentres()[facei]
                << " is the first occurrence."
                << exit(FatalError);
        }
    }
}


// To combineReduce a labelList. Filters out duplicates.
class uniqueEqOp
{

public:

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
    const UIndirectList<label>& extrudeMeshFaces,// mesh faces that are extruded
    const label meshEdgeI                       // mesh edge
)
{
    // Make set of extruded faces.
    labelHashSet extrudeFaceSet(extrudeMeshFaces.size());
    forAll(extrudeMeshFaces, i)
    {
        extrudeFaceSet.insert(extrudeMeshFaces[i]);
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
    const UIndirectList<label>& extrudeMeshFaces,// mesh faces that are extruded
    const label meshEdgeI                       // mesh edge
)
{
    // Make set of extruded faces.
    labelHashSet extrudeFaceSet(extrudeMeshFaces.size());
    forAll(extrudeMeshFaces, i)
    {
        extrudeFaceSet.insert(extrudeMeshFaces[i]);
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


// Calculate per edge min and max zone
void calcEdgeMinMaxZone
(
    const fvMesh& mesh,
    const primitiveFacePatch& extrudePatch,
    const labelList& extrudeMeshEdges,
    const labelList& zoneID,
    const mapDistribute& extrudeEdgeFacesMap,
    const labelListList& extrudeEdgeGlobalFaces,

    labelList& minZoneID,
    labelList& maxZoneID
)
{
    // Get zoneIDs in extrudeEdgeGlobalFaces order
    labelList mappedZoneID(zoneID);
    extrudeEdgeFacesMap.distribute(mappedZoneID);

    // Get min and max zone per edge
    minZoneID.setSize(extrudeEdgeGlobalFaces.size(), labelMax);
    maxZoneID.setSize(extrudeEdgeGlobalFaces.size(), labelMin);

    forAll(extrudeEdgeGlobalFaces, edgeI)
    {
        const labelList& eFaces = extrudeEdgeGlobalFaces[edgeI];
        if (eFaces.size())
        {
            forAll(eFaces, i)
            {
                label zoneI = mappedZoneID[eFaces[i]];
                minZoneID[edgeI] = min(minZoneID[edgeI], zoneI);
                maxZoneID[edgeI] = max(maxZoneID[edgeI], zoneI);
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh,
        extrudeMeshEdges,
        minZoneID,
        minEqOp<label>(),
        labelMax        // null value
    );
    syncTools::syncEdgeList
    (
        mesh,
        extrudeMeshEdges,
        maxZoneID,
        maxEqOp<label>(),
        labelMin        // null value
    );
}


// Count the number of faces in patches that need to be created. Calculates:
//  zoneSidePatch[zoneI]         : the number of side faces to be created
//  zoneZonePatch[zoneA,zoneB]   : the number of faces in between zoneA and B
// Since this only counts we're not taking the processor patches into
// account.
void countExtrudePatches
(
    const fvMesh& mesh,
    const label nZones,
    const primitiveFacePatch& extrudePatch,
    const labelList& extrudeMeshFaces,
    const labelList& extrudeMeshEdges,

    const labelListList& extrudeEdgeGlobalFaces,
    const labelList& minZoneID,
    const labelList& maxZoneID,

    labelList& zoneSidePatch,
    labelList& zoneZonePatch
)
{
    // Check on master edge for use of zones. Since we only want to know
    // whether they are being used at all no need to accurately count on slave
    // edge as well. Just add all together at the end of this routine so it
    // gets detected at least.

    forAll(extrudePatch.edgeFaces(), edgeI)
    {
        const labelList& eFaces = extrudePatch.edgeFaces()[edgeI];

        if (eFaces.size() == 2)
        {
            // Internal edge - check if in between different zones.
            if (minZoneID[edgeI] != maxZoneID[edgeI])
            {
                zoneZonePatch[minZoneID[edgeI]*nZones+maxZoneID[edgeI]]++;
            }
        }
        else if
        (
            eFaces.size() == 1
         && extrudeEdgeGlobalFaces[edgeI].size() == 2
        )
        {
            // Coupled edge - check if in between different zones.
            if (minZoneID[edgeI] != maxZoneID[edgeI])
            {
                const edge& e = extrudePatch.edges()[edgeI];
                const pointField& pts = extrudePatch.localPoints();
                WarningInFunction
                    << "Edge " << edgeI
                    << "at " << pts[e[0]] << pts[e[1]]
                    << " is a coupled edge and in between two different zones "
                    << minZoneID[edgeI] << " and " << maxZoneID[edgeI] << endl
                    << "    This is currently not supported." << endl;

                zoneZonePatch[minZoneID[edgeI]*nZones+maxZoneID[edgeI]]++;
            }
        }
        else
        {
            // One or more than two edge-faces.
            // Check whether we are on a mesh edge with external patches. If
            // so choose any uncovered one. If none found put face in
            // undetermined zone 'side' patch

            label facei = findUncoveredPatchFace
            (
                mesh,
                UIndirectList<label>(extrudeMeshFaces, eFaces),
                extrudeMeshEdges[edgeI]
            );

            if (facei == -1)
            {
                zoneSidePatch[minZoneID[edgeI]]++;
            }
        }
    }
    // Synchronise decistion. Actual numbers are not important, just make
    // sure that they're > 0 on all processors.
    Pstream::listCombineGather(zoneSidePatch, plusEqOp<label>());
    Pstream::listCombineScatter(zoneSidePatch);
    Pstream::listCombineGather(zoneZonePatch, plusEqOp<label>());
    Pstream::listCombineScatter(zoneZonePatch);
}


void addCouplingPatches
(
    const fvMesh& mesh,
    const word& regionName,
    const word& shellRegionName,
    const wordList& zoneNames,
    const wordList& zoneShadowNames,
    const boolList& isInternal,
    const labelList& zoneIDs,

    DynamicList<polyPatch*>& newPatches,
    labelList& interRegionTopPatch,
    labelList& interRegionBottomPatch
)
{
    Pout<< "Adding coupling patches:" << nl << nl
        << "patchID\tpatch\ttype" << nl
        << "-------\t-----\t----"
        << endl;

    interRegionTopPatch.setSize(zoneNames.size(), -1);
    interRegionBottomPatch.setSize(zoneNames.size(), -1);

    label nOldPatches = newPatches.size();
    forAll(zoneNames, zoneI)
    {
        word interName
        (
            regionName
           +"_to_"
           +shellRegionName
           +'_'
           +zoneNames[zoneI]
        );

        if (isInternal[zoneI])
        {
            interRegionTopPatch[zoneI] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName + "_top",
                newPatches
            );
            Pout<< interRegionTopPatch[zoneI]
                << '\t' << newPatches[interRegionTopPatch[zoneI]]->name()
                << '\t' << newPatches[interRegionTopPatch[zoneI]]->type()
                << nl;

            interRegionBottomPatch[zoneI] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName + "_bottom",
                newPatches
            );
            Pout<< interRegionBottomPatch[zoneI]
                << '\t' << newPatches[interRegionBottomPatch[zoneI]]->name()
                << '\t' << newPatches[interRegionBottomPatch[zoneI]]->type()
                << nl;
        }
        else if (zoneShadowNames.size() == 0)
        {
            interRegionTopPatch[zoneI] = addPatch<polyPatch>
            (
                mesh.boundaryMesh(),
                zoneNames[zoneI] + "_top",
                newPatches
            );
            Pout<< interRegionTopPatch[zoneI]
                << '\t' << newPatches[interRegionTopPatch[zoneI]]->name()
                << '\t' << newPatches[interRegionTopPatch[zoneI]]->type()
                << nl;

            interRegionBottomPatch[zoneI] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName,
                newPatches
            );
            Pout<< interRegionBottomPatch[zoneI]
                << '\t' << newPatches[interRegionBottomPatch[zoneI]]->name()
                << '\t' << newPatches[interRegionBottomPatch[zoneI]]->type()
                << nl;
        }
        else    // patch using shadow face zones.
        {
            interRegionTopPatch[zoneI] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                zoneShadowNames[zoneI] + "_top",
                newPatches
            );
            Pout<< interRegionTopPatch[zoneI]
                << '\t' << newPatches[interRegionTopPatch[zoneI]]->name()
                << '\t' << newPatches[interRegionTopPatch[zoneI]]->type()
                << nl;

            interRegionBottomPatch[zoneI] = addPatch<mappedWallPolyPatch>
            (
                mesh.boundaryMesh(),
                interName,
                newPatches
            );
            Pout<< interRegionBottomPatch[zoneI]
                << '\t' << newPatches[interRegionBottomPatch[zoneI]]->name()
                << '\t' << newPatches[interRegionBottomPatch[zoneI]]->type()
                << nl;
        }
    }
    Pout<< "Added " << newPatches.size()-nOldPatches
        << " inter-region patches." << nl
        << endl;
}


// Sets sidePatch[edgeI] to interprocessor or cyclic patch. Adds any
// coupled patches if necessary.
void addCoupledPatches
(
    const fvMesh& mesh,
    const primitiveFacePatch& extrudePatch,
    const labelList& extrudeMeshFaces,
    const labelList& extrudeMeshEdges,
    const mapDistribute& extrudeEdgeFacesMap,
    const labelListList& extrudeEdgeGlobalFaces,

    labelList& sidePatchID,
    DynamicList<polyPatch*>& newPatches
)
{
    // Calculate opposite processor for coupled edges (only if shared by
    // two procs). Note: could have saved original globalEdgeFaces structure.

    // Get procID in extrudeEdgeGlobalFaces order
    labelList procID(extrudeEdgeGlobalFaces.size(), Pstream::myProcNo());
    extrudeEdgeFacesMap.distribute(procID);

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
        extrudeMeshEdges,
        minProcID,
        minEqOp<label>(),
        labelMax        // null value
    );
    syncTools::syncEdgeList
    (
        mesh,
        extrudeMeshEdges,
        maxProcID,
        maxEqOp<label>(),
        labelMin        // null value
    );

    Pout<< "Adding processor or cyclic patches:" << nl << nl
        << "patchID\tpatch" << nl
        << "-------\t-----"
        << endl;

    label nOldPatches = newPatches.size();

    sidePatchID.setSize(extrudePatch.edgeFaces().size(), -1);
    forAll(extrudePatch.edgeFaces(), edgeI)
    {
        const labelList& eFaces = extrudePatch.edgeFaces()[edgeI];

        if
        (
            eFaces.size() == 1
         && extrudeEdgeGlobalFaces[edgeI].size() == 2
        )
        {
            // coupled boundary edge. Find matching patch.
            label nbrProci = minProcID[edgeI];
            if (nbrProci == Pstream::myProcNo())
            {
                nbrProci = maxProcID[edgeI];
            }


            if (nbrProci == Pstream::myProcNo())
            {
                // Cyclic patch since both procs the same. This cyclic should
                // already exist in newPatches so no adding necessary.

                label facei = findUncoveredCyclicPatchFace
                (
                    mesh,
                    UIndirectList<label>(extrudeMeshFaces, eFaces),
                    extrudeMeshEdges[edgeI]
                );

                if (facei != -1)
                {
                    const polyBoundaryMesh& patches = mesh.boundaryMesh();

                    label newPatchi = findPatchID
                    (
                        newPatches,
                        patches[patches.whichPatch(facei)].name()
                    );

                    sidePatchID[edgeI] = newPatchi;
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
                word name
                (
                    processorPolyPatch::newName(Pstream::myProcNo(), nbrProci)
                );

                sidePatchID[edgeI] = findPatchID(newPatches, name);

                if (sidePatchID[edgeI] == -1)
                {
                    dictionary patchDict;
                    patchDict.add("myProcNo", Pstream::myProcNo());
                    patchDict.add("neighbProcNo", nbrProci);

                    sidePatchID[edgeI] = addPatch<processorPolyPatch>
                    (
                        mesh.boundaryMesh(),
                        name,
                        patchDict,
                        newPatches
                    );

                    Pout<< sidePatchID[edgeI] << '\t' << name
                        << nl;
                }
            }
        }
    }
    Pout<< "Added " << newPatches.size()-nOldPatches
        << " coupled patches." << nl
        << endl;
}


void addZoneSidePatches
(
    const fvMesh& mesh,
    const wordList& zoneNames,
    const word& oneDPolyPatchType,

    DynamicList<polyPatch*>& newPatches,
    labelList& zoneSidePatch
)
{
    Pout<< "Adding patches for sides on zones:" << nl << nl
        << "patchID\tpatch" << nl
        << "-------\t-----"
        << endl;

    label nOldPatches = newPatches.size();

    forAll(zoneSidePatch, zoneI)
    {
        if (oneDPolyPatchType != word::null)
        {
            // Reuse single empty patch.
            word patchName;
            if (oneDPolyPatchType == "empty")
            {
                patchName = "oneDEmptyPatch";
                zoneSidePatch[zoneI] = addPatch<emptyPolyPatch>
                (
                    mesh.boundaryMesh(),
                    patchName,
                    newPatches
                );
            }
            else if (oneDPolyPatchType == "wedge")
            {
                patchName = "oneDWedgePatch";
                zoneSidePatch[zoneI] = addPatch<wedgePolyPatch>
                (
                    mesh.boundaryMesh(),
                    patchName,
                    newPatches
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Type " << oneDPolyPatchType << " does not exist "
                    << exit(FatalError);
            }

            Pout<< zoneSidePatch[zoneI] << '\t' << patchName << nl;
        }
        else if (zoneSidePatch[zoneI] > 0)
        {
            word patchName = zoneNames[zoneI] + "_" + "side";

            zoneSidePatch[zoneI] = addPatch<polyPatch>
            (
                mesh.boundaryMesh(),
                patchName,
                newPatches
            );

            Pout<< zoneSidePatch[zoneI] << '\t' << patchName << nl;
        }
    }
    Pout<< "Added " << newPatches.size()-nOldPatches << " zone-side patches."
        << nl << endl;
}


void addInterZonePatches
(
    const fvMesh& mesh,
    const wordList& zoneNames,
    const bool oneD,

    labelList& zoneZonePatch_min,
    labelList& zoneZonePatch_max,
    DynamicList<polyPatch*>& newPatches
)
{
    Pout<< "Adding inter-zone patches:" << nl << nl
        << "patchID\tpatch" << nl
        << "-------\t-----"
        << endl;

    dictionary transformDict;
    transformDict.add
    (
        "transform",
        cyclicPolyPatch::transformTypeNames[cyclicPolyPatch::NOORDERING]
    );

    label nOldPatches = newPatches.size();

    if (!oneD)
    {
        forAll(zoneZonePatch_min, minZone)
        {
            for (label maxZone = minZone; maxZone < zoneNames.size(); maxZone++)
            {
                label index = minZone*zoneNames.size()+maxZone;

                if (zoneZonePatch_min[index] > 0)
                {
                    word minToMax =
                        zoneNames[minZone]
                      + "_to_"
                      + zoneNames[maxZone];
                    word maxToMin =
                        zoneNames[maxZone]
                      + "_to_"
                      + zoneNames[minZone];

                    {
                        transformDict.set("neighbourPatch", maxToMin);
                        zoneZonePatch_min[index] =
                        addPatch<nonuniformTransformCyclicPolyPatch>
                        (
                            mesh.boundaryMesh(),
                            minToMax,
                            transformDict,
                            newPatches
                        );
                        Pout<< zoneZonePatch_min[index] << '\t' << minToMax
                            << nl;
                    }
                    {
                        transformDict.set("neighbourPatch", minToMax);
                        zoneZonePatch_max[index] =
                        addPatch<nonuniformTransformCyclicPolyPatch>
                        (
                            mesh.boundaryMesh(),
                            maxToMin,
                            transformDict,
                            newPatches
                        );
                        Pout<< zoneZonePatch_max[index] << '\t' << maxToMin
                            << nl;
                    }

                }
            }
        }
    }
    Pout<< "Added " << newPatches.size()-nOldPatches << " inter-zone patches."
        << nl << endl;
}


tmp<pointField> calcOffset
(
    const primitiveFacePatch& extrudePatch,
    const createShellMesh& extruder,
    const polyPatch& pp
)
{
    vectorField::subField fc = pp.faceCentres();

    tmp<pointField> toffsets(new pointField(fc.size()));
    pointField& offsets = toffsets.ref();

    forAll(fc, i)
    {
        label meshFacei = pp.start()+i;
        label patchFacei = mag(extruder.faceToFaceMap()[meshFacei])-1;
        point patchFc = extrudePatch[patchFacei].centre
        (
            extrudePatch.points()
        );
        offsets[i] = patchFc - fc[i];
    }
    return toffsets;
}


void setCouplingInfo
(
    fvMesh& mesh,
    const labelList& zoneToPatch,
    const word& sampleRegion,
    const mappedWallPolyPatch::sampleMode mode,
    const List<pointField>& offsets
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    List<polyPatch*> newPatches
    (
        patches.size(),
        static_cast<polyPatch*>(nullptr)
    );

    forAll(zoneToPatch, zoneI)
    {
        label patchi = zoneToPatch[zoneI];

        if (patchi != -1)
        {
            const polyPatch& pp = patches[patchi];

            if (isA<mappedWallPolyPatch>(pp))
            {
                newPatches[patchi] = new mappedWallPolyPatch
                (
                    pp.name(),
                    pp.size(),
                    pp.start(),
                    patchi,
                    sampleRegion,                           // sampleRegion
                    mode,                                   // sampleMode
                    pp.name(),                              // samplePatch
                    offsets[zoneI],                         // offset
                    patches
                );
            }
        }
    }

    forAll(newPatches, patchi)
    {
        if (!newPatches[patchi])
        {
            newPatches[patchi] = patches[patchi].clone(patches).ptr();
        }
    }

    mesh.removeFvBoundary();
    mesh.addFvPatches(newPatches, true);
}


// Extrude and write geometric properties
void extrudeGeometricProperties
(
    const polyMesh& mesh,
    const primitiveFacePatch& extrudePatch,
    const createShellMesh& extruder,
    const polyMesh& regionMesh,
    const extrudeModel& model
)
{
     const pointIOField patchFaceCentres
     (
        IOobject
        (
            "patchFaceCentres",
            mesh.pointsInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );

    const pointIOField patchEdgeCentres
    (
        IOobject
        (
            "patchEdgeCentres",
            mesh.pointsInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        )
    );

    // forAll(extrudePatch.edges(), edgeI)
    //{
    //    const edge& e = extrudePatch.edges()[edgeI];
    //    Pout<< "Edge:" << e.centre(extrudePatch.localPoints()) << nl
    //        << "read:" << patchEdgeCentres[edgeI]
    //        << endl;
    //}


    // Determine edge normals on original patch
    labelList patchEdges;
    labelList coupledEdges;
    PackedBoolList sameEdgeOrientation;
    PatchTools::matchEdges
    (
        extrudePatch,
        mesh.globalData().coupledPatch(),
        patchEdges,
        coupledEdges,
        sameEdgeOrientation
    );

    pointField patchEdgeNormals
    (
        PatchTools::edgeNormals
        (
            mesh,
            extrudePatch,
            patchEdges,
            coupledEdges
        )
    );


    pointIOField faceCentres
    (
        IOobject
        (
            "faceCentres",
            regionMesh.pointsInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        regionMesh.nFaces()
    );


    // Work out layers. Guaranteed in columns so no fancy parallel bits.


    forAll(extruder.faceToFaceMap(), facei)
    {
        if (extruder.faceToFaceMap()[facei] != 0)
        {
            // 'horizontal' face
            label patchFacei = mag(extruder.faceToFaceMap()[facei])-1;

            label celli = regionMesh.faceOwner()[facei];
            if (regionMesh.isInternalFace(facei))
            {
                celli = max(celli, regionMesh.faceNeighbour()[facei]);
            }

            // Calculate layer from cell numbering (see createShellMesh)
            label layerI = (celli % model.nLayers());

            if
            (
               !regionMesh.isInternalFace(facei)
             && extruder.faceToFaceMap()[facei] > 0
            )
            {
                // Top face
                layerI++;
            }


            // Recalculate based on extrusion model
            faceCentres[facei] = model
            (
                patchFaceCentres[patchFacei],
                extrudePatch.faceNormals()[patchFacei],
                layerI
            );
        }
        else
        {
            // 'vertical face
            label patchEdgeI = extruder.faceToEdgeMap()[facei];
            label layerI =
            (
                regionMesh.faceOwner()[facei]
              % model.nLayers()
            );

            // Extrude patch edge centre to this layer
            point pt0 = model
            (
                patchEdgeCentres[patchEdgeI],
                patchEdgeNormals[patchEdgeI],
                layerI
            );
            // Extrude patch edge centre to next layer
            point pt1 = model
            (
                patchEdgeCentres[patchEdgeI],
                patchEdgeNormals[patchEdgeI],
                layerI+1
            );

            // Interpolate
            faceCentres[facei] = 0.5*(pt0+pt1);
        }
    }

    pointIOField cellCentres
    (
        IOobject
        (
            "cellCentres",
            regionMesh.pointsInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        regionMesh.nCells()
    );

    forAll(extruder.cellToFaceMap(), celli)
    {
        label patchFacei = extruder.cellToFaceMap()[celli];

        // Calculate layer from cell numbering (see createShellMesh)
        label layerI = (celli % model.nLayers());

        // Recalculate based on extrusion model
        point pt0 = model
        (
            patchFaceCentres[patchFacei],
            extrudePatch.faceNormals()[patchFacei],
            layerI
        );
        point pt1 = model
        (
            patchFaceCentres[patchFacei],
            extrudePatch.faceNormals()[patchFacei],
            layerI+1
        );

        // Interpolate
        cellCentres[celli] = 0.5*(pt0+pt1);
    }


    // Bit of checking
    if (false)
    {
        OBJstream faceStr(regionMesh.time().path()/"faceCentres.obj");
        OBJstream cellStr(regionMesh.time().path()/"cellCentres.obj");

        forAll(faceCentres, facei)
        {
            Pout<< "Model     :" << faceCentres[facei] << endl
                << "regionMesh:" << regionMesh.faceCentres()[facei] << endl;
            faceStr.write
            (
                linePointRef
                (
                    faceCentres[facei],
                    regionMesh.faceCentres()[facei]
                )
            );
        }
        forAll(cellCentres, celli)
        {
            Pout<< "Model     :" << cellCentres[celli] << endl
                << "regionMesh:" << regionMesh.cellCentres()[celli] << endl;
            cellStr.write
            (
                linePointRef
                (
                    cellCentres[celli],
                    regionMesh.cellCentres()[celli]
                )
            );
        }
    }



    Info<< "Writing geometric properties for mesh " << regionMesh.name()
        << " to " << regionMesh.pointsInstance() << nl
        << endl;

    bool ok = faceCentres.write() && cellCentres.write();

    if (!ok)
    {
        FatalErrorInFunction
            << "Failed writing " << faceCentres.objectPath()
            << " and " << cellCentres.objectPath()
            << exit(FatalError);
    }
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


    const word oldInstance = mesh.pointsInstance();
    bool overwrite = args.optionFound("overwrite");


    const word dictName("extrudeToRegionMeshDict");

    #include "setSystemMeshDictionaryIO.H"

    IOdictionary dict(dictIO);


    // Point generator
    autoPtr<extrudeModel> model(extrudeModel::New(dict));


    // Region
    const word shellRegionName(dict.lookup("region"));

    // Faces to extrude - either faceZones or faceSets (boundary faces only)
    wordList zoneNames;
    wordList zoneShadowNames;

    bool hasZones = dict.found("faceZones");
    if (hasZones)
    {
        dict.lookup("faceZones") >> zoneNames;
        dict.readIfPresent("faceZonesShadow", zoneShadowNames);

        // Check
        if (dict.found("faceSets"))
        {
            FatalIOErrorIn(args.executable().c_str(), dict)
                << "Please supply faces to extrude either through 'faceZones'"
                << " or 'faceSets' entry. Found both."
                << exit(FatalIOError);
        }
    }
    else
    {
        dict.lookup("faceSets") >> zoneNames;
        dict.readIfPresent("faceSetsShadow", zoneShadowNames);
    }


    mappedPatchBase::sampleMode sampleMode =
        mappedPatchBase::sampleModeNames_[dict.lookup("sampleMode")];

    const Switch oneD(dict.lookup("oneD"));
    Switch oneDNonManifoldEdges(false);
    word oneDPatchType(emptyPolyPatch::typeName);
    if (oneD)
    {
        oneDNonManifoldEdges = dict.lookupOrDefault("nonManifold", false);
        dict.lookup("oneDPolyPatchType") >> oneDPatchType;
    }

    const Switch adaptMesh(dict.lookup("adaptMesh"));

    if (hasZones)
    {
        Info<< "Extruding zones " << zoneNames
            << " on mesh " << regionName
            << " into shell mesh " << shellRegionName
            << endl;
    }
    else
    {
        Info<< "Extruding faceSets " << zoneNames
            << " on mesh " << regionName
            << " into shell mesh " << shellRegionName
            << endl;
    }

    if (shellRegionName == regionName)
    {
        FatalIOErrorIn(args.executable().c_str(), dict)
            << "Cannot extrude into same region as mesh." << endl
            << "Mesh region : " << regionName << endl
            << "Shell region : " << shellRegionName
            << exit(FatalIOError);
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




    //// Read objects in time directory
    // IOobjectList objects(mesh, runTime.timeName());
    //
    //// Read vol fields.
    //
    // PtrList<volScalarField> vsFlds;
    // ReadFields(mesh, objects, vsFlds);
    //
    // PtrList<volVectorField> vvFlds;
    // ReadFields(mesh, objects, vvFlds);
    //
    // PtrList<volSphericalTensorField> vstFlds;
    // ReadFields(mesh, objects, vstFlds);
    //
    // PtrList<volSymmTensorField> vsymtFlds;
    // ReadFields(mesh, objects, vsymtFlds);
    //
    // PtrList<volTensorField> vtFlds;
    // ReadFields(mesh, objects, vtFlds);
    //
    //// Read surface fields.
    //
    // PtrList<surfaceScalarField> ssFlds;
    // ReadFields(mesh, objects, ssFlds);
    //
    // PtrList<surfaceVectorField> svFlds;
    // ReadFields(mesh, objects, svFlds);
    //
    // PtrList<surfaceSphericalTensorField> sstFlds;
    // ReadFields(mesh, objects, sstFlds);
    //
    // PtrList<surfaceSymmTensorField> ssymtFlds;
    // ReadFields(mesh, objects, ssymtFlds);
    //
    // PtrList<surfaceTensorField> stFlds;
    // ReadFields(mesh, objects, stFlds);
    //
    //// Read point fields.
    //
    // PtrList<pointScalarField> psFlds;
    // ReadFields(pointMesh::New(mesh), objects, psFlds);
    //
    // PtrList<pointVectorField> pvFlds;
    // ReadFields(pointMesh::New(mesh), objects, pvFlds);



    // Create dummy fv* files
    createDummyFvMeshFiles(mesh, shellRegionName);


    word meshInstance;
    if (!overwrite)
    {
        runTime++;
        meshInstance = runTime.timeName();
    }
    else
    {
        meshInstance = oldInstance;
    }
    Info<< "Writing meshes to " << meshInstance << nl << endl;


    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Extract faces to extrude
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: zoneID are regions of extrusion. They are not mesh.faceZones
    // indices.

    // From extrude zone to mesh zone (or -1 if extruding faceSets)
    labelList meshZoneID;
    // Per extrude zone whether contains internal or external faces
    boolList isInternal(zoneNames.size(), false);

    labelList extrudeMeshFaces;
    faceList zoneFaces;
    labelList zoneID;
    boolList zoneFlipMap;
    // Shadow
    labelList zoneShadowIDs;    // from extrude shadow zone to mesh zone
    labelList extrudeMeshShadowFaces;
    boolList zoneShadowFlipMap;
    labelList zoneShadowID;

    if (hasZones)
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        meshZoneID.setSize(zoneNames.size());
        forAll(zoneNames, i)
        {
            meshZoneID[i] = faceZones.findZoneID(zoneNames[i]);
            if (meshZoneID[i] == -1)
            {
                FatalIOErrorIn(args.executable().c_str(), dict)
                    << "Cannot find zone " << zoneNames[i] << endl
                    << "Valid zones are " << faceZones.names()
                    << exit(FatalIOError);
            }
        }
        // Collect per face information
        label nExtrudeFaces = 0;
        forAll(meshZoneID, i)
        {
            nExtrudeFaces += faceZones[meshZoneID[i]].size();
        }
        extrudeMeshFaces.setSize(nExtrudeFaces);
        zoneFaces.setSize(nExtrudeFaces);
        zoneID.setSize(nExtrudeFaces);
        zoneFlipMap.setSize(nExtrudeFaces);
        nExtrudeFaces = 0;
        forAll(meshZoneID, i)
        {
            const faceZone& fz = faceZones[meshZoneID[i]];
            const primitiveFacePatch& fzp = fz();
            forAll(fz, j)
            {
                extrudeMeshFaces[nExtrudeFaces] = fz[j];
                zoneFaces[nExtrudeFaces] = fzp[j];
                zoneID[nExtrudeFaces] = i;
                zoneFlipMap[nExtrudeFaces] = fz.flipMap()[j];
                nExtrudeFaces++;

                if (mesh.isInternalFace(fz[j]))
                {
                    isInternal[i] = true;
                }
            }
        }

        // Shadow zone
        // ~~~~~~~~~~~

        if (zoneShadowNames.size())
        {
            zoneShadowIDs.setSize(zoneShadowNames.size());
            forAll(zoneShadowNames, i)
            {
                zoneShadowIDs[i] = faceZones.findZoneID(zoneShadowNames[i]);
                if (zoneShadowIDs[i] == -1)
                {
                    FatalIOErrorIn(args.executable().c_str(), dict)
                        << "Cannot find zone " << zoneShadowNames[i] << endl
                        << "Valid zones are " << faceZones.names()
                        << exit(FatalIOError);
                }
            }

            label nShadowFaces = 0;
            forAll(zoneShadowIDs, i)
            {
                nShadowFaces += faceZones[zoneShadowIDs[i]].size();
            }

            extrudeMeshShadowFaces.setSize(nShadowFaces);
            zoneShadowFlipMap.setSize(nShadowFaces);
            zoneShadowID.setSize(nShadowFaces);

            nShadowFaces = 0;
            forAll(zoneShadowIDs, i)
            {
                const faceZone& fz = faceZones[zoneShadowIDs[i]];
                forAll(fz, j)
                {
                    extrudeMeshShadowFaces[nShadowFaces] = fz[j];
                    zoneShadowFlipMap[nShadowFaces] = fz.flipMap()[j];
                    zoneShadowID[nShadowFaces] = i;
                    nShadowFaces++;
                }
            }
        }
    }
    else
    {
        meshZoneID.setSize(zoneNames.size(), -1);
        // Load faceSets
        PtrList<faceSet> zones(zoneNames.size());
        forAll(zoneNames, i)
        {
            Info<< "Loading faceSet " << zoneNames[i] << endl;
            zones.set(i, new faceSet(mesh, zoneNames[i]));
        }


        // Collect per face information
        label nExtrudeFaces = 0;
        forAll(zones, i)
        {
            nExtrudeFaces += zones[i].size();
        }
        extrudeMeshFaces.setSize(nExtrudeFaces);
        zoneFaces.setSize(nExtrudeFaces);
        zoneID.setSize(nExtrudeFaces);
        zoneFlipMap.setSize(nExtrudeFaces);

        nExtrudeFaces = 0;
        forAll(zones, i)
        {
            const faceSet& fz = zones[i];
            forAllConstIter(faceSet, fz, iter)
            {
                label facei = iter.key();
                if (mesh.isInternalFace(facei))
                {
                    FatalIOErrorIn(args.executable().c_str(), dict)
                        << "faceSet " << fz.name()
                        << "contains internal faces."
                        << " This is not permitted."
                        << exit(FatalIOError);
                }
                extrudeMeshFaces[nExtrudeFaces] = facei;
                zoneFaces[nExtrudeFaces] = mesh.faces()[facei];
                zoneID[nExtrudeFaces] = i;
                zoneFlipMap[nExtrudeFaces] = false;
                nExtrudeFaces++;

                if (mesh.isInternalFace(facei))
                {
                    isInternal[i] = true;
                }
            }
        }


        // Shadow zone
        // ~~~~~~~~~~~

        PtrList<faceSet> shadowZones(zoneShadowNames.size());
        if (zoneShadowNames.size())
        {
            zoneShadowIDs.setSize(zoneShadowNames.size(), -1);
            forAll(zoneShadowNames, i)
            {
                shadowZones.set(i, new faceSet(mesh, zoneShadowNames[i]));
            }

            label nShadowFaces = 0;
            forAll(shadowZones, i)
            {
                nShadowFaces += shadowZones[i].size();
            }

            if (nExtrudeFaces != nShadowFaces)
            {
                FatalIOErrorIn(args.executable().c_str(), dict)
                    << "Extruded faces " << nExtrudeFaces << endl
                    << "is different from shadow faces. " << nShadowFaces
                    << "This is not permitted " << endl
                    << exit(FatalIOError);
            }

            extrudeMeshShadowFaces.setSize(nShadowFaces);
            zoneShadowFlipMap.setSize(nShadowFaces);
            zoneShadowID.setSize(nShadowFaces);

            nShadowFaces = 0;
            forAll(shadowZones, i)
            {
                const faceSet& fz = shadowZones[i];
                forAllConstIter(faceSet, fz, iter)
                {
                    label facei = iter.key();
                    if (mesh.isInternalFace(facei))
                    {
                        FatalIOErrorIn(args.executable().c_str(), dict)
                            << "faceSet " << fz.name()
                            << "contains internal faces."
                            << " This is not permitted."
                            << exit(FatalIOError);
                    }
                    extrudeMeshShadowFaces[nShadowFaces] = facei;
                    zoneShadowFlipMap[nShadowFaces] = false;
                    zoneShadowID[nShadowFaces] = i;
                    nShadowFaces++;
                }
            }
        }
    }
    const primitiveFacePatch extrudePatch(zoneFaces.xfer(), mesh.points());


    Pstream::listCombineGather(isInternal, orEqOp<bool>());
    Pstream::listCombineScatter(isInternal);

    // Check zone either all internal or all external faces
    checkZoneInside(mesh, zoneNames, zoneID, extrudeMeshFaces, isInternal);



    const pointField& extrudePoints = extrudePatch.localPoints();
    const faceList& extrudeFaces = extrudePatch.localFaces();
    const labelListList& edgeFaces = extrudePatch.edgeFaces();


    Info<< "extrudePatch :"
        << " faces:" << extrudePatch.size()
        << " points:" << extrudePatch.nPoints()
        << " edges:" << extrudePatch.nEdges()
        << nl
        << endl;


    // Determine per-extrude-edge info
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Corresponding mesh edges
    const labelList extrudeMeshEdges
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
            extrudeMeshEdges
        )
    );
    List<Map<label>> compactMap;
    const mapDistribute extrudeEdgeFacesMap
    (
        globalExtrudeFaces,
        extrudeEdgeGlobalFaces,
        compactMap
    );


    // Determine min and max zone per edge
    labelList edgeMinZoneID;
    labelList edgeMaxZoneID;
    calcEdgeMinMaxZone
    (
        mesh,
        extrudePatch,
        extrudeMeshEdges,
        zoneID,
        extrudeEdgeFacesMap,
        extrudeEdgeGlobalFaces,

        edgeMinZoneID,
        edgeMaxZoneID
    );




    DynamicList<polyPatch*> regionPatches(patches.size());
    // Copy all non-local patches since these are used on boundary edges of
    // the extrusion
    forAll(patches, patchi)
    {
        if (!isA<processorPolyPatch>(patches[patchi]))
        {
            label newPatchi = regionPatches.size();
            regionPatches.append
            (
                patches[patchi].clone
                (
                    patches,
                    newPatchi,
                    0,              // size
                    0               // start
                ).ptr()
            );
        }
    }


    // Add interface patches
    // ~~~~~~~~~~~~~~~~~~~~~

    // From zone to interface patch (region side)
    labelList interRegionTopPatch;
    labelList interRegionBottomPatch;

    addCouplingPatches
    (
        mesh,
        regionName,
        shellRegionName,
        zoneNames,
        zoneShadowNames,
        isInternal,
        meshZoneID,

        regionPatches,
        interRegionTopPatch,
        interRegionBottomPatch
    );


    // From zone to interface patch (mesh side)
    labelList interMeshTopPatch;
    labelList interMeshBottomPatch;

    if (adaptMesh)
    {
        // Add coupling patches to mesh

        // Clone existing patches
        DynamicList<polyPatch*> newPatches(patches.size());
        forAll(patches, patchi)
        {
            newPatches.append(patches[patchi].clone(patches).ptr());
        }

        // Add new patches
        addCouplingPatches
        (
            mesh,
            regionName,
            shellRegionName,
            zoneNames,
            zoneShadowNames,
            isInternal,
            meshZoneID,

            newPatches,
            interMeshTopPatch,
            interMeshBottomPatch
        );

        // Add to mesh
        mesh.clearOut();
        mesh.removeFvBoundary();
        mesh.addFvPatches(newPatches, true);

        //!Note: from this point on mesh patches differs from regionPatches
    }


    // Patch per extruded face
    labelList extrudeTopPatchID(extrudePatch.size());
    labelList extrudeBottomPatchID(extrudePatch.size());

    forAll(zoneID, facei)
    {
        extrudeTopPatchID[facei] = interRegionTopPatch[zoneID[facei]];
        extrudeBottomPatchID[facei] = interRegionBottomPatch[zoneID[facei]];
    }



    // Count how many patches on special edges of extrudePatch are necessary
    // - zoneXXX_sides
    // - zoneXXX_zoneYYY
    labelList zoneSidePatch(zoneNames.size(), 0);
    // Patch to use for minZone
    labelList zoneZonePatch_min(zoneNames.size()*zoneNames.size(), 0);
    // Patch to use for maxZone
    labelList zoneZonePatch_max(zoneNames.size()*zoneNames.size(), 0);

    countExtrudePatches
    (
        mesh,
        zoneNames.size(),

        extrudePatch,           // patch
        extrudeMeshFaces,       // mesh face per patch face
        extrudeMeshEdges,       // mesh edge per patch edge

        extrudeEdgeGlobalFaces, // global indexing per patch edge
        edgeMinZoneID,          // minZone per patch edge
        edgeMaxZoneID,          // maxZone per patch edge

        zoneSidePatch,          // per zone-side num edges that extrude into it
        zoneZonePatch_min       // per zone-zone num edges that extrude into it
    );

    // Now we'll have:
    //  zoneSidePatch[zoneA] : number of faces needed on the side of zoneA
    //  zoneZonePatch_min[zoneA,zoneB] : number of faces needed in between A,B


    // Add the zone-side patches.
    addZoneSidePatches
    (
        mesh,
        zoneNames,
        (oneD ? oneDPatchType : word::null),

        regionPatches,
        zoneSidePatch
    );


    // Add the patches in between zones
    addInterZonePatches
    (
        mesh,
        zoneNames,
        oneD,

        zoneZonePatch_min,
        zoneZonePatch_max,
        regionPatches
    );


    // Sets sidePatchID[edgeI] to interprocessor patch. Adds any
    // interprocessor or cyclic patches if necessary.
    labelList sidePatchID;
    addCoupledPatches
    (
        mesh,
        extrudePatch,
        extrudeMeshFaces,
        extrudeMeshEdges,
        extrudeEdgeFacesMap,
        extrudeEdgeGlobalFaces,

        sidePatchID,
        regionPatches
    );


//    // Add all the newPatches to the mesh and fields
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    {
//        forAll(newPatches, patchi)
//        {
//            Pout<< "Adding patch " << patchi
//                << " name:" << newPatches[patchi]->name()
//                << endl;
//        }
//        // label nOldPatches = mesh.boundary().size();
//        mesh.clearOut();
//        mesh.removeFvBoundary();
//        mesh.addFvPatches(newPatches, true);
//        //// Add calculated fvPatchFields for the added patches
//        // for
//        //(
//        //    label patchi = nOldPatches;
//        //    patchi < mesh.boundary().size();
//        //    patchi++
//        //)
//        //{
//        //    Pout<< "ADDing calculated to patch " << patchi
//        //        << endl;
//        //    addCalculatedPatchFields(mesh);
//        //}
//        // Pout<< "** Added " << mesh.boundary().size()-nOldPatches
//        //    << " patches." << endl;
//    }


    // Set patches to use for edges to be extruded into boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // In order of edgeFaces: per edge, per originating face the
    // patch to use for the side face (from the extruded edge).
    // If empty size create an internal face.
    labelListList extrudeEdgePatches(extrudePatch.nEdges());

    // Is edge a non-manifold edge
    PackedBoolList nonManifoldEdge(extrudePatch.nEdges());

    // Note: logic has to be same as in countExtrudePatches.
    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        labelList& ePatches = extrudeEdgePatches[edgeI];

        if (oneD)
        {
            ePatches.setSize(eFaces.size());
            forAll(eFaces, i)
            {
                ePatches[i] = zoneSidePatch[zoneID[eFaces[i]]];
            }

            if (oneDNonManifoldEdges)
            {
                //- Set nonManifoldEdge[edgeI] for non-manifold edges only
                //  The other option is to have non-manifold edges everywhere
                //  and generate space overlapping columns of cells.
                if (eFaces.size() != 2)
                {
                    nonManifoldEdge[edgeI] = 1;
                }
            }
            else
            {
                nonManifoldEdge[edgeI] = 1;
            }
        }
        else if (eFaces.size() == 2)
        {
            label zone0 = zoneID[eFaces[0]];
            label zone1 = zoneID[eFaces[1]];

            if (zone0 != zone1) // || (cos(angle) > blabla))
            {
                label minZone = min(zone0,zone1);
                label maxZone = max(zone0,zone1);
                label index = minZone*zoneNames.size()+maxZone;

                ePatches.setSize(eFaces.size());

                if (zone0 == minZone)
                {
                    ePatches[0] = zoneZonePatch_min[index];
                    ePatches[1] = zoneZonePatch_max[index];
                }
                else
                {
                    ePatches[0] = zoneZonePatch_max[index];
                    ePatches[1] = zoneZonePatch_min[index];
                }

                nonManifoldEdge[edgeI] = 1;
            }
        }
        else if (sidePatchID[edgeI] != -1)
        {
            // Coupled extrusion
            ePatches.setSize(eFaces.size());
            forAll(eFaces, i)
            {
                ePatches[i] = sidePatchID[edgeI];
            }
        }
        else
        {
            label facei = findUncoveredPatchFace
            (
                mesh,
                UIndirectList<label>(extrudeMeshFaces, eFaces),
                extrudeMeshEdges[edgeI]
            );

            if (facei != -1)
            {
                label newPatchi = findPatchID
                (
                    regionPatches,
                    patches[patches.whichPatch(facei)].name()
                );
                ePatches.setSize(eFaces.size(), newPatchi);
            }
            else
            {
                ePatches.setSize(eFaces.size());
                forAll(eFaces, i)
                {
                    ePatches[i] = zoneSidePatch[zoneID[eFaces[i]]];
                }
            }
            nonManifoldEdge[edgeI] = 1;
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
                label localRegionI = pRegions[fp];
                localSum[localRegionI] += extrudePatch.faceNormals()[facei];
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
        localRegionNormals /= mag(localRegionNormals);
    }


    // For debugging: dump hedgehog plot of normals
    if (false)
    {
        OFstream str(runTime.path()/"localRegionNormals.obj");
        label vertI = 0;

        scalar thickness = model().sumThickness(1);

        forAll(pointLocalRegions, facei)
        {
            const face& f = extrudeFaces[facei];

            forAll(f, fp)
            {
                label region = pointLocalRegions[facei][fp];
                const point& pt = extrudePoints[f[fp]];

                meshTools::writeOBJ(str, pt);
                vertI++;
                meshTools::writeOBJ
                (
                    str,
                    pt+thickness*localRegionNormals[region]
                );
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Use model to create displacements of first layer
    vectorField firstDisp(localRegionNormals.size());
    forAll(firstDisp, regionI)
    {
        // const point& regionPt = regionCentres[regionI];
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

    createShellMesh extruder
    (
        extrudePatch,
        pointLocalRegions,
        localRegionPoints
    );


    autoPtr<mapPolyMesh> shellMap;
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
        xferCopy(pointField()),
        xferCopy(faceList()),
        xferCopy(labelList()),
        xferCopy(labelList()),
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

    {
        polyTopoChange meshMod(regionPatches.size());

        extruder.setRefinement
        (
            firstDisp,                              // first displacement
            model().expansionRatio(),
            model().nLayers(),                      // nLayers
            extrudeTopPatchID,
            extrudeBottomPatchID,
            extrudeEdgePatches,
            meshMod
        );

        // Enforce actual point posititions according to extrudeModel (model)
        // (extruder.setRefinement only does fixed expansionRatio)
        // The regionPoints and nLayers are looped in the same way as in
        // createShellMesh
        DynamicList<point>& newPoints = const_cast<DynamicList<point>&>
        (
            meshMod.points()
        );
        label meshPointi = extrudePatch.localPoints().size();
        forAll(localRegionPoints, regionI)
        {
            label pointi = localRegionPoints[regionI];
            point pt = extrudePatch.localPoints()[pointi];
            const vector& n = localRegionNormals[regionI];

            for (label layerI = 1; layerI <= model().nLayers(); layerI++)
            {
                newPoints[meshPointi++] = model()(pt, n, layerI);
            }
        }

        shellMap = meshMod.changeMesh
        (
            regionMesh,     // mesh to change
            false           // inflate
        );
    }

    // Necessary?
    regionMesh.setInstance(meshInstance);


    // Update numbering on extruder.
    extruder.updateMesh(shellMap);


    // Calculate offsets from shell mesh back to original mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<pointField> topOffsets(zoneNames.size());
    List<pointField> bottomOffsets(zoneNames.size());

    forAll(regionMesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = regionMesh.boundaryMesh()[patchi];

        if (isA<mappedWallPolyPatch>(pp))
        {
            if (findIndex(interRegionTopPatch, patchi) != -1)
            {
                label zoneI = findIndex(interRegionTopPatch, patchi);
                topOffsets[zoneI] = calcOffset(extrudePatch, extruder, pp);
            }
            else if (findIndex(interRegionBottomPatch, patchi) != -1)
            {
                label zoneI = findIndex(interRegionBottomPatch, patchi);
                bottomOffsets[zoneI] = calcOffset(extrudePatch, extruder, pp);
            }
        }
    }


    // Change top and bottom boundary conditions on regionMesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Correct top patches for offset
        setCouplingInfo
        (
            regionMesh,
            interRegionTopPatch,
            regionName,                 // name of main mesh
            sampleMode,                 // sampleMode
            topOffsets
        );

        // Correct bottom patches for offset
        setCouplingInfo
        (
            regionMesh,
            interRegionBottomPatch,
            regionName,
            sampleMode,                 // sampleMode
            bottomOffsets
        );

        // Remove any unused patches
        deleteEmptyPatches(regionMesh);
    }

    // Change top and bottom boundary conditions on main mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (adaptMesh)
    {
        // Correct top patches for offset
        setCouplingInfo
        (
            mesh,
            interMeshTopPatch,
            shellRegionName,                        // name of shell mesh
            sampleMode,                             // sampleMode
            -topOffsets
        );

        // Correct bottom patches for offset
        setCouplingInfo
        (
            mesh,
            interMeshBottomPatch,
            shellRegionName,
            sampleMode,
            -bottomOffsets
        );
    }



    // Write addressing from region mesh back to originating patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelIOList cellToPatchFaceAddressing
    (
        IOobject
        (
            "cellToPatchFaceAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.cellToFaceMap()
    );
    cellToPatchFaceAddressing.note() = "cell to patch face addressing";

    labelIOList faceToPatchFaceAddressing
    (
        IOobject
        (
            "faceToPatchFaceAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.faceToFaceMap()
    );
    faceToPatchFaceAddressing.note() =
        "front/back face + turning index to patch face addressing";

    labelIOList faceToPatchEdgeAddressing
    (
        IOobject
        (
            "faceToPatchEdgeAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.faceToEdgeMap()
    );
    faceToPatchEdgeAddressing.note() =
        "side face to patch edge addressing";

    labelIOList pointToPatchPointAddressing
    (
        IOobject
        (
            "pointToPatchPointAddressing",
            regionMesh.facesInstance(),
            regionMesh.meshSubDir,
            regionMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        extruder.pointToPointMap()
    );
    pointToPatchPointAddressing.note() =
        "point to patch point addressing";


    Info<< "Writing mesh " << regionMesh.name()
        << " to " << regionMesh.facesInstance() << nl
        << endl;

    bool ok =
        regionMesh.write()
     && cellToPatchFaceAddressing.write()
     && faceToPatchFaceAddressing.write()
     && faceToPatchEdgeAddressing.write()
     && pointToPatchPointAddressing.write();

    if (!ok)
    {
        FatalErrorInFunction
            << "Failed writing mesh " << regionMesh.name()
            << " at location " << regionMesh.facesInstance()
            << exit(FatalError);
    }


    // See if we need to extrude coordinates as well
    {
        autoPtr<pointIOField> patchFaceCentresPtr;

        IOobject io
        (
            "patchFaceCentres",
            mesh.pointsInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::MUST_READ
        );
        if (io.typeHeaderOk<pointIOField>(true))
        {
            // Read patchFaceCentres and patchEdgeCentres
            Info<< "Reading patch face,edge centres : "
                << io.name() << " and patchEdgeCentres" << endl;

            extrudeGeometricProperties
            (
                mesh,
                extrudePatch,
                extruder,
                regionMesh,
                model()
            );
        }
    }




    // Insert baffles into original mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<mapPolyMesh> addBafflesMap;

    if (adaptMesh)
    {
        polyTopoChange meshMod(mesh);

        // Modify faces to be in bottom (= always coupled) patch
        forAll(extrudeMeshFaces, zoneFacei)
        {
            label meshFacei = extrudeMeshFaces[zoneFacei];
            label zoneI = zoneID[zoneFacei];
            bool flip = zoneFlipMap[zoneFacei];
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
                    interMeshBottomPatch[zoneI],// patch for face
                    meshZoneID[zoneI],          // zone for face
                    flip                        // face flip in zone
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
                    interMeshBottomPatch[zoneI],    // patch for face
                    meshZoneID[zoneI],              // zone for face
                    !flip                           // face flip in zone
                );
            }
        }

        if (zoneShadowNames.size() > 0) // if there is a top faceZone specified
        {
            forAll(extrudeMeshFaces, zoneFacei)
            {
                label meshFacei = extrudeMeshShadowFaces[zoneFacei];
                label zoneI = zoneShadowID[zoneFacei];
                bool flip = zoneShadowFlipMap[zoneFacei];
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
                        interMeshTopPatch[zoneI],   // patch for face
                        meshZoneID[zoneI],          // zone for face
                        flip                        // face flip in zone
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
                        interMeshTopPatch[zoneI],       // patch for face
                        meshZoneID[zoneI],              // zone for face
                        !flip                           // face flip in zone
                    );
                }
            }
        }
        else
        {
            // Add faces (using same points) to be in top patch
            forAll(extrudeMeshFaces, zoneFacei)
            {
                label meshFacei = extrudeMeshFaces[zoneFacei];
                label zoneI = zoneID[zoneFacei];
                bool flip = zoneFlipMap[zoneFacei];
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
                            interMeshTopPatch[zoneI],       // patch for face
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
                        interMeshTopPatch[zoneI],       // patch for face
                        -1,                             // zone for face
                        false                           // zone flip
                    );
                }
            }
        }

        // Change the mesh. Change points directly (no inflation).
        addBafflesMap = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(addBafflesMap);


//XXXXXX
// Update maps! e.g. faceToPatchFaceAddressing
//XXXXXX

        // Move mesh (since morphing might not do this)
        if (addBafflesMap().hasMotionPoints())
        {
            mesh.movePoints(addBafflesMap().preMotionPoints());
        }

        mesh.setInstance(meshInstance);

        // Remove any unused patches
        deleteEmptyPatches(mesh);

        Info<< "Writing mesh " << mesh.name()
            << " to " << mesh.facesInstance() << nl
            << endl;

        if (!mesh.write())
        {
            FatalErrorInFunction
                << "Failed writing mesh " << mesh.name()
                << " at location " << mesh.facesInstance()
                << exit(FatalError);
        }
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
