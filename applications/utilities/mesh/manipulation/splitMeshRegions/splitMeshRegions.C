/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    splitMeshRegions

Description
    Splits mesh into multiple region meshes.

Usage
    \b splitMeshRegions -disconnected -largest
    \b splitMeshRegions -cellZones all
    \b splitMeshRegions -cellZones 'element heatsink' -defaultRegion air

    Options:
      - \par -disconnected \n
        Put disconnected parts of the mesh into separate regions

      - \par -disconnectedFaceSet <name> \n
        Specify a set of faces that are considered disconnected

      - \par -cellZones <names> \n
        Put cells in the specified zones into separate regions

      - \par -defaultRegion <name> \n
        Name given to regions which do not correspond to a named cell zone.
        Defaults to "region". Multiple such regions will have an index appended
        (i.e., "region1", "region2", ...).

      - \par -largest \n
        Only write the largest region

      - \par -insidePoint <point> \n
        Only write the region containing the given point

      - \par -noMeshes \n
        Do not write region meshes

      - \par -noFields \n
        Do not write fields

      - \par -minOverlapFraction <fraction> \n
        The minimum fraction a zone must overlap a cell-zone in
        order to be named after it and associated with it. Defaults to 0.

      - \par -faceZonePatches \n
        Use face zones to define the inter-region patches instead of creating a
        single inter-region patch per pair of regions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"
#include "IOobjectList.H"
#include "faceSet.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "syncTools.H"
#include "ReadFields.H"
#include "mappedWallPolyPatch.H"
#include "fvMeshTools.H"
#include "meshSearch.H"
#include "RemoteData.H"
#include "IOmanip.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Standard default region base name
const word standardRegionName("region");


// Prepend prefix to selected patches.
void renamePatches
(
    fvMesh& mesh,
    const word& prefix,
    const labelList& patchesToRename
)
{
    polyBoundaryMesh& patches =
        const_cast<polyBoundaryMesh&>(mesh.poly().boundary());

    // Create a list of all the new names
    wordList newNames = patches.names();
    forAll(patchesToRename, i)
    {
        const label patchi = patchesToRename[i];
        newNames[patchi] = prefix + '_' + newNames[patchi];
    }

    // Apply to the patches
    patches.renamePatches(newNames, true);
}


template<class GeoField>
void subsetVolFields
(
    const fvMesh& mesh,
    const fvMesh& subMesh,
    const labelList& cellMap,
    const labelList& faceMap,
    const labelHashSet& addedPatches
)
{
    const labelList patchMap(identityMap(mesh.poly().boundary().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< indent << "Mapping " << GeoField::typeName
            << " " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                cellMap,
                faceMap
            )
        );

        // Initialise the introduced patch fields with the adjacent cell value
        forAll(tSubFld().boundaryField(), patchi)
        {
            if (addedPatches.found(patchi))
            {
                tSubFld.ref().boundaryFieldRef()[patchi] ==
                    tSubFld.ref().boundaryField()[patchi].patchInternalField();
            }
        }

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt() = IOobject::AUTO_WRITE;
        subFld->store();
    }
}


template<class GeoField>
void subsetSurfaceFields
(
    const fvMesh& mesh,
    const fvMesh& subMesh,
    const labelList& cellMap,
    const labelList& faceMap,
    const labelHashSet& addedPatches
)
{
    const labelList patchMap(identityMap(mesh.poly().boundary().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< indent << "Mapping " << GeoField::typeName
            << " " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                cellMap,
                faceMap
            )
        );

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt() = IOobject::AUTO_WRITE;
        subFld->store();
    }
}


// Select all cells not in the region
labelList getNonRegionCells(const labelList& cellRegion, const label regioni)
{
    DynamicList<label> nonRegionCells(cellRegion.size());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] != regioni)
        {
            nonRegionCells.append(celli);
        }
    }
    return nonRegionCells.shrink();
}


void addToInterface
(
    const polyMesh& mesh,
    const label zoneID,
    const label ownRegion,
    const label neiRegion,
    EdgeMap<Map<label>>& regionsToSize
)
{
    edge interface
    (
        min(ownRegion, neiRegion),
        max(ownRegion, neiRegion)
    );

    EdgeMap<Map<label>>::iterator iter = regionsToSize.find
    (
        interface
    );

    if (iter != regionsToSize.end())
    {
        // Check if zone present
        Map<label>::iterator zoneFnd = iter().find(zoneID);
        if (zoneFnd != iter().end())
        {
            // Found zone. Increment count.
            zoneFnd()++;
        }
        else
        {
            // New or no zone. Insert with count 1.
            iter().insert(zoneID, 1);
        }
    }
    else
    {
        // Create new interface of size 1.
        Map<label> zoneToSize;
        zoneToSize.insert(zoneID, 1);
        regionsToSize.insert(interface, zoneToSize);
    }
}


label whichZone
(
    const polyMesh& mesh,
    const bool faceZonePatches,
    const label facei
)
{
    if (faceZonePatches)
    {
        const labelList zones(mesh.faceZones().whichZones(facei));

        if (zones.size() == 0)
        {
            return -1;
        }
        else if (zones.size() == 1)
        {
            return zones[0];
        }
        else
        {
            FatalErrorInFunction
                << "Face " << facei << " is in more than one zone " << zones
                << exit(FatalError);

            return -1;
        }
    }
    else
    {
        return -1;
    }
}


// Get region-region interface name and sizes.
// Returns interfaces as straight list for looping in identical order.
void getInterfaceSizes
(
    const polyMesh& mesh,
    const bool faceZonePatches,
    const labelList& cellRegion,
    const wordList& regionNames,

    edgeList& interfaces,
    List<Pair<word>>& interfaceNames,
    labelList& interfaceNFaces,
    labelList& faceInterfaces
)
{
    // From region-region to faceZone (or -1) to number of faces.

    EdgeMap<Map<label>> regionsToSize;


    // Internal faces
    // ~~~~~~~~~~~~~~

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                whichZone(mesh, faceZonePatches, facei),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }

    // Boundary faces
    // ~~~~~~~~~~~~~~

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nFaces()-mesh.nInternalFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);

    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                whichZone(mesh, faceZonePatches, facei),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }


    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive and add to my sizes
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);

                EdgeMap<Map<label>> slaveSizes(fromSlave);

                forAllConstIter(EdgeMap<Map<label>>, slaveSizes, slaveIter)
                {
                    EdgeMap<Map<label>>::iterator masterIter =
                        regionsToSize.find(slaveIter.key());

                    if (masterIter != regionsToSize.end())
                    {
                        // Same inter-region
                        const Map<label>& slaveInfo = slaveIter();
                        Map<label>& masterInfo = masterIter();

                        forAllConstIter(Map<label>, slaveInfo, iter)
                        {
                            label zoneID = iter.key();
                            label slaveSize = iter();

                            Map<label>::iterator zoneFnd = masterInfo.find
                            (
                                zoneID
                            );
                            if (zoneFnd != masterInfo.end())
                            {
                                zoneFnd() += slaveSize;
                            }
                            else
                            {
                                masterInfo.insert(zoneID, slaveSize);
                            }
                        }
                    }
                    else
                    {
                        regionsToSize.insert(slaveIter.key(), slaveIter());
                    }
                }
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::blocking,
                    Pstream::masterNo()
                );
                toMaster << regionsToSize;
            }
        }
    }

    // Rework

    Pstream::scatter(regionsToSize);



    // Now we have the global sizes of all inter-regions.
    // Invert this on master and distribute.
    label nInterfaces = 0;
    forAllConstIter(EdgeMap<Map<label>>, regionsToSize, iter)
    {
        const Map<label>& info = iter();
        nInterfaces += info.size();
    }

    interfaces.setSize(nInterfaces);
    interfaceNames.setSize(nInterfaces);
    interfaceNFaces.setSize(nInterfaces);
    EdgeMap<Map<label>> regionsToInterface(nInterfaces);

    nInterfaces = 0;
    forAllConstIter(EdgeMap<Map<label>>, regionsToSize, iter)
    {
        const edge& e = iter.key();
        const word& name0 = regionNames[e[0]];
        const word& name1 = regionNames[e[1]];

        const Map<label>& info = iter();
        forAllConstIter(Map<label>, info, infoIter)
        {
            interfaces[nInterfaces] = iter.key();
            label zoneID = infoIter.key();
            if (zoneID == -1)
            {
                interfaceNames[nInterfaces] = Pair<word>
                (
                    name0 + "_to_" + name1,
                    name1 + "_to_" + name0
                );
            }
            else
            {
                const word& zoneName = mesh.faceZones()[zoneID].name();
                interfaceNames[nInterfaces] = Pair<word>
                (
                    zoneName + "_" + name0 + "_to_" + name1,
                    zoneName + "_" + name1 + "_to_" + name0
                );
            }
            interfaceNFaces[nInterfaces] = infoIter();

            if (regionsToInterface.found(e))
            {
                regionsToInterface[e].insert(zoneID, nInterfaces);
            }
            else
            {
                Map<label> zoneAndInterface;
                zoneAndInterface.insert(zoneID, nInterfaces);
                regionsToInterface.insert(e, zoneAndInterface);
            }
            nInterfaces++;
        }
    }


    // Now all processor have consistent interface information

    Pstream::scatter(interfaces);
    Pstream::scatter(interfaceNames);
    Pstream::scatter(interfaceNFaces);
    Pstream::scatter(regionsToInterface);

    // Mark all inter-region faces.
    faceInterfaces.setSize(mesh.nFaces(), -1);

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            const label zoneID = whichZone(mesh, faceZonePatches, facei);

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceInterfaces[facei] = regionsToInterface[interface][zoneID];
        }
    }
    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            const label zoneID = whichZone(mesh, faceZonePatches, facei);

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceInterfaces[facei] = regionsToInterface[interface][zoneID];
        }
    }
}


// Create mesh for region.
autoPtr<polyTopoChangeMap> createRegionMesh
(
    const fvMesh& mesh,
    // Region info
    const labelList& cellRegion,
    const label regioni,
    const word& regionName,
    // Interface info
    const labelList& interfacePatches,
    const labelList& faceInterfaces,

    autoPtr<fvMesh>& newMesh
)
{
    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nFaces()-mesh.nInternalFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);


    // Topology change container. Start off from existing mesh.
    polyTopoChange meshMod(mesh);

    // Cell remover engine
    removeCells cellRemover(mesh);

    // Select all but region cells
    labelList cellsToRemove(getNonRegionCells(cellRegion, regioni));

    // Find out which faces will get exposed. Note that this
    // gets faces in mesh face order. So both regions will get same
    // face in same order (important!)
    labelList exposedFaces = cellRemover.getExposedFaces(cellsToRemove);

    labelList exposedPatchIDs(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        label interfacei = faceInterfaces[facei];

        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = -1;

        if (mesh.isInternalFace(facei))
        {
            neiRegion = cellRegion[mesh.faceNeighbour()[facei]];
        }
        else
        {
            neiRegion = coupledRegion[facei-mesh.nInternalFaces()];
        }


        // Check which side is being kept - determines which of the two
        // patches will be used.

        label otherRegion = -1;

        if (ownRegion == regioni && neiRegion != regioni)
        {
            otherRegion = neiRegion;
        }
        else if (ownRegion != regioni && neiRegion == regioni)
        {
            otherRegion = ownRegion;
        }
        else
        {
            FatalErrorInFunction
                << "Exposed face:" << facei
                << " fc:" << mesh.faceCentres()[facei]
                << " has owner region " << ownRegion
                << " and neighbour region " << neiRegion
                << " when handling region:" << regioni
                << exit(FatalError);
        }

        // Find the patch.
        if (regioni < otherRegion)
        {
            exposedPatchIDs[i] = interfacePatches[interfacei];
        }
        else
        {
            exposedPatchIDs[i] = interfacePatches[interfacei]+1;
        }
    }

    // Remove faces
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    autoPtr<polyTopoChangeMap> map = meshMod.makeMesh
    (
        newMesh,
        IOobject
        (
            regionName,
            mesh.time().name(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    return map;
}


void createAndWriteRegion
(
    const fvMesh& mesh,
    const labelList& cellRegion,
    const wordList& regionNames,
    const labelList& faceInterfaces,
    const labelList& interfacePatches,
    const label regioni,
    const word& newMeshInstance
)
{
    Info<< "Creating mesh for region \'"
        << regionNames[regioni] << "\':" << endl << incrIndent;

    autoPtr<fvMesh> newMesh;
    autoPtr<polyTopoChangeMap> map = createRegionMesh
    (
        mesh,
        cellRegion,
        regioni,
        regionNames[regioni],
        interfacePatches,
        faceInterfaces,
        newMesh
    );

    // Make map of all added patches
    labelHashSet addedPatches(2*interfacePatches.size());
    forAll(interfacePatches, interfacei)
    {
        addedPatches.insert(interfacePatches[interfacei]);
        addedPatches.insert(interfacePatches[interfacei]+1);
    }

    Info<< indent << "Mapping fields" << endl << incrIndent;

    // Map existing fields
    newMesh().topoChange(map());

    // Add subsetted fields
    #define SUBSET_VOL_FIELDS(Type, nullArg) \
        subsetVolFields<VolField<Type>>      \
        (                                    \
            mesh,                            \
            newMesh(),                       \
            map().cellMap(),                 \
            map().faceMap(),                 \
            addedPatches                     \
        );
    FOR_ALL_FIELD_TYPES(SUBSET_VOL_FIELDS);
    #undef SUBSET_VOL_FIELDS

    #define SUBSET_SURFACE_FIELDS(Type, nullArg) \
        subsetSurfaceFields<SurfaceField<Type>>  \
        (                                        \
            mesh,                                \
            newMesh(),                           \
            map().cellMap(),                     \
            map().faceMap(),                     \
            addedPatches                         \
        );
    FOR_ALL_FIELD_TYPES(SUBSET_SURFACE_FIELDS);
    #undef SUBSET_SURFACE_FIELDS

    Info<< decrIndent;

    const polyBoundaryMesh& newPatches = newMesh().poly().boundary();
    newPatches.checkParallelSync(true);

    // Delete empty patches ...

    // Create reordering list to move patches-to-be-deleted to end
    labelList oldToNew(newPatches.size(), -1);
    DynamicList<label> sharedPatches(newPatches.size());
    label newI = 0;

    Info<< indent << "Deleting empty patches" << endl;

    // Assumes all non-proc boundaries are on all processors!
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if (returnReduce(pp.size(), sumOp<label>()) > 0)
            {
                oldToNew[patchi] = newI;
                if (!addedPatches.found(patchi))
                {
                    sharedPatches.append(newI);
                }
                newI++;
            }
        }
    }

    // Same for processor patches (but need no reduction)
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (isA<processorPolyPatch>(pp) && pp.size())
        {
            oldToNew[patchi] = newI++;
        }
    }

    const label nNewPatches = newI;

    // Move all delete-able patches to the end
    forAll(oldToNew, patchi)
    {
        if (oldToNew[patchi] == -1)
        {
            oldToNew[patchi] = newI++;
        }
    }

    fvMeshTools::reorderPatches(newMesh(), oldToNew, nNewPatches, true);

    Info<< indent << "Writing new mesh" << endl;

    newMesh().setInstance(newMeshInstance);
    newMesh().write();

    Info<< endl << decrIndent;
}


// Create for every region-region interface with non-zero size two patches.
// First one is for minimumregion to maximumregion.
// Note that patches get created in same order on all processors (if parallel)
// since looping over synchronised 'interfaces'.
labelList addRegionPatches
(
    fvMesh& mesh,
    const wordList& regionNames,
    const edgeList& interfaces,
    const List<Pair<word>>& interfaceNames
)
{
    labelList interfacePatches(interfaces.size());

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];
        const Pair<word>& names = interfaceNames[interI];

        mappedWallPolyPatch patch1
        (
            names[0],
            0,                  // overridden
            0,                  // overridden
            0,                  // overridden
            regionNames[e[1]],  // neighbourRegion
            names[1],           // neighbourPatch
            mesh.poly().boundary()
        );

        interfacePatches[interI] = fvMeshTools::addPatch(mesh, patch1);

        mappedWallPolyPatch patch2
        (
            names[1],
            0,
            0,
            0,
            regionNames[e[0]],  // neighbourRegion
            names[0],
            mesh.poly().boundary()
        );

        fvMeshTools::addPatch(mesh, patch2);
    }

    fvMeshTools::addedPatches(mesh);

    return interfacePatches;
}


void getCellCellZoneis
(
    const polyMesh& mesh,
    const labelList cellZoneis,
    labelList& cellCellZones,
    labelList& bFaceNbrCellZones
)
{
    cellCellZones.setSize(mesh.nCells());
    cellCellZones = -1;

    label failCelli = -1;
    DynamicList<label> failCellZones;

    forAll(cellZoneis, i)
    {
        const label cellZonei = cellZoneis[i];

        const cellZone& cz = mesh.cellZones()[cellZonei];

        forAll(cz, czCelli)
        {
            const label celli = cz[czCelli];

            if (cellCellZones[celli] == -1)
            {
                cellCellZones[celli] = cellZonei;
            }
            else
            {
                if (failCelli == -1)
                {
                    failCelli = celli;
                    failCellZones.append(cellCellZones[celli]);
                }
                if (failCelli == celli)
                {
                    failCellZones.append(cellZonei);
                }
            }
        }
    }

    if (failCelli != -1)
    {
        FatalErrorInFunction
            << "Cell " << failCelli << " with centre "
            << mesh.cellCentres()[failCelli]
            << " is multiple selected zones; ";
        forAll(failCellZones, i)
        {
            FatalError<< mesh.cellZones()[failCellZones[i]].name();
            if (i < failCellZones.size() - 2) FatalError << ", ";
            if (i == failCellZones.size() - 2) FatalError << " and ";
        }
        FatalError
            << ". This is not allowed."
            << exit(FatalError);
    }

    bFaceNbrCellZones.setSize(mesh.nFaces() - mesh.nInternalFaces());
    bFaceNbrCellZones = -1;

    forAll(bFaceNbrCellZones, bFacei)
    {
        bFaceNbrCellZones[bFacei] =
            cellCellZones[mesh.faceOwner()[mesh.nInternalFaces() + bFacei]];
    }

    syncTools::swapBoundaryFaceList(mesh, bFaceNbrCellZones);
}


void matchRegions
(
    const polyMesh& mesh,
    const label nRegions,
    const labelList& cellRegion,
    const word& defaultRegion,
    labelList& cellZoneRegions,
    wordList& regionNames,
    labelList& regionCellZones,
    const scalar minOverlapFraction
)
{
    cellZoneRegions.setSize(mesh.cellZones().size(), -1);
    regionNames.setSize(nRegions);
    regionCellZones.setSize(nRegions, -1);

    // Check the zone names are in sync
    {
        List<wordList> cellZoneNames(Pstream::nProcs());

        cellZoneNames[Pstream::myProcNo()] = mesh.cellZones().toc();

        Pstream::gatherList(cellZoneNames);
        Pstream::scatterList(cellZoneNames);

        forAll(cellZoneNames, proci)
        {
            if (cellZoneNames[proci] != cellZoneNames[0])
            {
                FatalErrorInFunction
                    << "cellZones not synchronised across processors." << endl
                    << "Master has cellZones " << cellZoneNames[0] << endl
                    << "Processor " << proci << " has cellZones "
                    << cellZoneNames[proci] << exit(FatalError);
            }
        }
    }

    // Get the global numbers of cells in the zones
    labelList cellZoneNCells(mesh.cellZones().size(), 0);
    {
        forAll(mesh.cellZones(), zoneI)
        {
            cellZoneNCells[zoneI] = returnReduce
            (
                mesh.cellZones()[zoneI].size(),
                sumOp<label>()
            );
        }
    }

    boolList cellInZone(mesh.nCells(), false);

    Info<< "Trying to match regions to existing cell zones:" << endl;

    forAll(mesh.cellZones(), cellZonei)
    {
        UIndirectList<bool> cellZoneCellInZone
        (
            cellInZone,
            mesh.cellZones()[cellZonei]
        );

        cellZoneCellInZone = true;

        const label minOverlapNCells =
            max(1, label(minOverlapFraction*cellZoneNCells[cellZonei]));

        // Per region the number of cells in the zone
        labelList nCellsInZone(nRegions, 0);
        forAll(cellRegion, celli)
        {
            if (cellInZone[celli])
            {
                nCellsInZone[cellRegion[celli]]++;
            }
        }
        Pstream::listCombineGather(nCellsInZone, plusEqOp<label>());
        Pstream::listCombineScatter(nCellsInZone);

        // Pick the region with largest overlap of the zone
        label regioni = findMax(nCellsInZone);

        if (nCellsInZone[regioni] < minOverlapNCells)
        {
            // This region covers too little of zone. Not good enough.
            regioni = -1;
        }
        else
        {
            // Check that region contains no cells outside of the zone
            forAll(cellRegion, celli)
            {
                if (cellRegion[celli] == regioni && !cellInZone[celli])
                {
                    regioni = -1;
                    break;
                }
            }
            reduce(regioni, minOp<label>());
        }

        // If successful, name the identified region and map it to the zone
        if (regioni != -1)
        {
            Info<< "    Matched region " << regioni
                << " to zone " << mesh.cellZones()[cellZonei].name()
                << " with " << cellZoneNCells[cellZonei] << " cells "
                << endl;

            if (cellZoneRegions[cellZonei] == -1)
            {
                cellZoneRegions[cellZonei] = regioni;
                regionNames[regioni] = mesh.cellZones()[cellZonei].name();
                regionCellZones[regioni] = cellZonei;
            }
            else
            {
                cellZoneRegions[cellZonei] = -2;
                regionNames[regioni] = word::null;
                regionCellZones[regioni] = -1;
            }
        }

        cellZoneCellInZone = false;
    }

    // Restore the standard "this cell zone has no region" index
    forAll(cellZoneRegions, cellZonei)
    {
        cellZoneRegions[cellZonei] = max(cellZoneRegions[cellZonei], -1);
    }

    // Check if there are any regions which were not matched
    label nUnmatchedRegions = 0;
    forAll(regionCellZones, regioni)
    {
        if (regionCellZones[regioni] == -1)
        {
            nUnmatchedRegions ++;
        }
    }

    if (nUnmatchedRegions)
    {
        label nUnmatchedi = 1;

        // Create default names for unmatches regions
        forAll(regionCellZones, regioni)
        {
            if (regionCellZones[regioni] == -1)
            {
                if
                (
                    nUnmatchedRegions == 1
                 && defaultRegion != standardRegionName
                )
                {
                    regionNames[regioni] = defaultRegion;
                }
                else
                {
                    regionNames[regioni] =
                        defaultRegion + Foam::name(nUnmatchedi ++);
                }
            }
        }
    }

    Info<< endl;
}


void writeCellToRegion(const fvMesh& mesh, const labelList& cellRegion)
{
    // Write to manual decomposition option
    {
        labelIOList cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.facesInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            cellRegion
        );
        cellToRegion.write();

        Info<< "Writing region index per cell as a "
            << labelIOList::typeName << " to "
            << cellToRegion.relativeObjectPath() << nl << endl;
    }

    // Write for post processing
    {
        volScalarField::Internal cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, 0)
        );
        forAll(cellRegion, celli)
        {
            cellToRegion[celli] = cellRegion[celli];
        }
        cellToRegion.write();

        Info<< "Writing region index per cell as a "
            << volScalarField::Internal::typeName << " to "
            << cellToRegion.relativeObjectPath() << nl << endl;
    }
}


template<unsigned int NColumns>
void printTable(const List<FixedList<string, NColumns>>& table)
{
    FixedList<string::size_type, NColumns> tableWs(string::size_type(0));
    forAll(table, i)
    {
        forAll(table[i], j)
        {
            tableWs[j] = max(tableWs[j], table[i][j].size());
        }
    }

    forAll(table, i)
    {
        if (i == 1)
        {
            forAll(table[i], j)
            {
                if (j) Info << "  ";
                Info<< string(tableWs[j], '-').c_str();
            }
            Info<< endl;
        }

        forAll(table[i], j)
        {
            if (j) Info << "  ";
            Info<< setw(tableWs[j]) << table[i][j].c_str();
        }
        Info<< nl;
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "splits mesh into multiple regions (detected by walking across faces)"
    );
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    #include "addNoOverwriteOption.H"
    argList::addBoolOption
    (
        "disconnected",
        "put disconnected parts of the mesh into separate regions"
    );
    argList::addOption
    (
        "disconnectedFaceSet",
        "name",
        "specify a set of faces that are considered disconnected"
    );
    argList::addOption
    (
        "cellZones",
        "names",
        "put cells in the specified zones into separate regions"
    );
    argList::addOption
    (
        "defaultRegion",
        "name",
        "name given to regions which do not correspond to a named cell zone; "
        "default \"region\""
    );
    argList::addBoolOption
    (
        "largest",
        "only write the largest region"
    );
    argList::addOption
    (
        "insidePoint",
        "point",
        "only write the region containing the given point"
    );
    argList::addBoolOption
    (
        "noMeshes",
        "do not write region meshes"
    );
    argList::addBoolOption
    (
        "noFields",
        "do not write region fields"
    );
    argList::addOption
    (
        "minOverlapFraction",
        "fraction",
        "the minimum fraction a zone must overlap a cell-zone in order to "
        "be named after it and associated with it; default 0"
    );
    argList::addBoolOption
    (
        "faceZonePatches",
        "Use face zones to define the inter-region patches instead of creating "
        "a single inter-region patch per pair of regions"
    );

    #include "setRootCaseNoFunctionObjects.H"

    if (!args.optionFound("disconnected") && !args.optionFound("cellZones"))
    {
        FatalErrorInFunction
            << "Neither -disconnected nor -cellZones specified, so no "
            << "criteria is available to identify regions" << exit(FatalError);
    }

    #include "createTimeNoFunctionObjects.H"
    #include "createSpecifiedMeshNoChangers.H"

    const word oldInstance = mesh.pointsInstance();

    const bool largest = args.optionFound("largest");
    const bool insidePoint = args.optionFound("insidePoint");
    if (largest && insidePoint)
    {
        FatalErrorInFunction
            << "-largest (i.e., write the region with most cells)"
            << " and -insidePoint (write the region containing the given "
            << "point) cannot be used together" << exit(FatalError);
    }

    #include "setNoOverwrite.H"

    const bool faceZonePatches = args.optionFound("faceZonePatches");
    if (faceZonePatches)
    {
        Info<< "Using face zones to divide inter-region interfaces"
            << " into multiple patches" << nl << endl;
    }
    else
    {
        Info<< "Creating single patch per inter-region interface"
            << nl << endl;
    }

    const word defaultRegion
    (
        args.optionLookupOrDefault("defaultRegion", standardRegionName)
    );

    // Read the set of cell zones (if any) to operate on
    labelList cellZoneis;
    if (args.optionFound("cellZones"))
    {
        const wordList cellZoneNames =
            args.optionReadList<word>("cellZones");

        DynamicList<label> cellZoneisDyn(cellZoneNames.size());
        boolList cellZoneUsed(mesh.cellZones().size(), false);
        forAll(cellZoneNames, i)
        {
            if (cellZoneNames[i] == "all")
            {
                cellZoneisDyn = identityMap(mesh.cellZones().size());
                break;
            }

            const label cellZonei =
                mesh.cellZones().findIndex(cellZoneNames[i]);

            if (cellZonei == -1)
            {
                FatalErrorInFunction
                    << "Cell zone \'" << cellZoneNames[i] << "\' not found. "
                    << "Available cell zones are: " << mesh.cellZones().toc()
                    << exit(FatalError);
            }

            if (cellZoneUsed[cellZonei]) continue;

            cellZoneisDyn.append(cellZonei);
            cellZoneUsed[cellZonei] = true;
        }
        cellZoneis.transfer(cellZoneisDyn);
    }

    // Identify the regions. For every cell, assign a region index. Also, where
    // possible, create correspondence between the existing cellZones and the
    // new regions. Create region names.
    labelList cellRegions(mesh.nCells(), -1);
    labelList cellZoneRegions(mesh.cellZones().size(), -1);
    wordList regionNames;
    labelList regionCellZones;
    if (args.optionFound("disconnected"))
    {
        boolList faceBlockeds(mesh.nFaces(), false);

        const word disconnectedFaceSetName =
            args.optionLookupOrDefault("disconnectedFaceSet", word::null);

        if (disconnectedFaceSetName != word::null)
        {
            faceSet blockedFaceSet(mesh, disconnectedFaceSetName);

            Info<< "Read "
                << returnReduce(blockedFaceSet.size(), sumOp<label>())
                << " blocked faces from set \'" << disconnectedFaceSetName
                << "\'" << nl << endl;

            forAllConstIter(faceSet, blockedFaceSet, iter)
            {
                faceBlockeds[iter.key()] = true;
            }
        }

        if (args.optionFound("cellZones"))
        {
            labelList cellZoneRegions(mesh.cellZones().size(), -1);
            forAll(cellZoneis, i)
            {
                cellZoneRegions[cellZoneis[i]] = i;
            }

            labelList cellCellZones, bFaceNbrCellZones;
            getCellCellZoneis
            (
                mesh,
                cellZoneis,
                cellCellZones,
                bFaceNbrCellZones
            );

            forAll(mesh.faces(), facei)
            {
                const label ownerZone =
                    cellCellZones[mesh.faceOwner()[facei]];
                const label neighbourZone =
                    facei < mesh.nInternalFaces()
                  ? cellCellZones[mesh.faceNeighbour()[facei]]
                  : bFaceNbrCellZones[facei - mesh.nInternalFaces()];

                const label ownerRegion =
                    ownerZone != -1 ? cellZoneRegions[ownerZone] : -1;
                const label neighbourRegion =
                    neighbourZone != -1 ? cellZoneRegions[neighbourZone] : -1;

                faceBlockeds[facei] =
                    faceBlockeds[facei] || ownerRegion != neighbourRegion;
            }
        }

        regionSplit regions(mesh, faceBlockeds);

        cellRegions.transfer(regions);

        matchRegions
        (
            mesh,
            regions.nRegions(),
            cellRegions,
            defaultRegion,
            cellZoneRegions,
            regionNames,
            regionCellZones,
            args.optionLookupOrDefault<scalar>("minOverlapFraction", 0)
        );
    }
    else // args.optionFound("cellZones")
    {
        regionNames.resize(cellZoneis.size());
        regionCellZones.resize(cellZoneis.size());

        forAll(cellZoneis, i)
        {
            cellZoneRegions[cellZoneis[i]] = i;

            regionNames[i] = mesh.cellZones()[cellZoneis[i]].name();
            regionCellZones[i] = cellZoneis[i];
        }

        labelList cellCellZones, bFaceNbrCellZones;
        getCellCellZoneis
        (
            mesh,
            cellZoneis,
            cellCellZones,
            bFaceNbrCellZones
        );

        renumber(cellZoneRegions, cellCellZones, cellRegions);

        bool haveDefaultRegion = false;
        forAll(cellRegions, celli)
        {
            if (cellRegions[celli] == -1)
            {
                cellRegions[celli] = cellZoneis.size();
                haveDefaultRegion = true;
            }
        }
        reduce(haveDefaultRegion, orOp<bool>());

        if (haveDefaultRegion)
        {
            regionNames.append(defaultRegion);
            regionCellZones.append(-1);
        }
    }

    // Quit if there is only one region
    if (regionNames.size() == 1)
    {
        Info<< "Only one region identified. Doing nothing." << nl << endl;
        Info<< "End" << nl << endl;
        return 0;
    }

    // Calculate the number of cells in each region
    labelList regionNCells(regionNames.size(), 0);
    {
        forAll(cellRegions, celli)
        {
            regionNCells[cellRegions[celli]] ++;
        }

        forAll(regionNCells, regioni)
        {
            reduce(regionNCells[regioni], sumOp<label>());
        }
    }

    // Report identified regions
    {
        List<FixedList<string, 3>> table(regionNames.size() + 1);
        table[0] = {"Region", "nCells", "cellZone"};
        forAll(regionNCells, regioni)
        {
            table[regioni + 1][0] = regionNames[regioni];
            table[regioni + 1][1] = Foam::name(regionNCells[regioni]);
            table[regioni + 1][2] =
                regionCellZones[regioni] == -1
              ? "<none>"
              : mesh.cellZones()[regionCellZones[regioni]].name();
        }
        printTable(table);
        Info<< endl;
    }

    // Write decomposition to file
    writeCellToRegion(mesh, cellRegions);

    // Everything from here relates to writing the region meshes, so quit now
    // if we are not doing that
    if (args.optionFound("noMeshes"))
    {
        Info<< "End" << nl << endl;
        return 0;
    }

    // Make sure patches and zones are synchronised
    mesh.poly().boundary().checkParallelSync(true);
    mesh.faceZones().checkParallelSync(true);

    // Create information regarding interfaces. Each interface is identified by
    // the indices of the two regions on either side. This is stored as an edge
    // so that an EdgeMap can then be used later (...). Also constructed are
    // the names of the interfaces, the global number of faces they contain,
    // and a map from each mesh face to the corresponding interface (or -1).
    edgeList interfaces;
    List<Pair<word>> interfaceNames;
    labelList interfaceNFaces;
    labelList faceInterfaces;
    getInterfaceSizes
    (
        mesh,
        faceZonePatches,
        cellRegions,
        regionNames,
        interfaces,
        interfaceNames,
        interfaceNFaces,
        faceInterfaces
    );

    // Report interfaces
    if (interfaces.size())
    {
        List<FixedList<string, 4>> table(interfaces.size() + 1);
        table[0] =
        {
            "Interface",
            "Owner Region",
            "Neighbour Region",
            "nFaces"
        };
        forAll(interfaces, interfacei)
        {
            const edge& e = interfaces[interfacei];

            table[interfacei + 1][0] = interfaceNames[interfacei][0];
            table[interfacei + 1][1] = regionNames[e[0]];
            table[interfacei + 1][2] = regionNames[e[1]];
            table[interfacei + 1][3] = Foam::name(interfaceNFaces[interfacei]);
        }
        printTable(table);
        Info<< endl;
    }

    // Read objects in time directory
    const bool fields = !args.optionFound("noFields");
    IOobjectList objects(mesh, runTime.name());
    if (fields) Info<< "Reading geometric fields" << endl << incrIndent;
    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    // #include "readPointFields.H"
    if (fields) Info<< endl << decrIndent;

    // Remove any demand-driven fields (e.g., cell-centres)
    mesh.clearOut();

    // Add patches for interfaces. These get added to the current mesh, so that
    // they are then inherited by the region subsets. All combinations are
    // added, and the patches are initially empty. Patches will then be
    // populated within the region meshes, and empty patches will be filtered
    // out and removed afterwards.
    if (interfaces.size())
    {
        Info<< "Adding interface patches" << nl << endl;
    }
    else
    {
        Info<< "No interfaces to add" << nl << endl;
    }
    const labelList interfacePatches
    (
        addRegionPatches
        (
            mesh,
            regionNames,
            interfaces,
            interfaceNames
        )
    );

    if (!overwrite)
    {
        runTime++;
    }

    // Create regions
    if (largest)
    {
        const label regioni = findMax(regionNCells);

        Info<< "Subsetting largest region " << regioni << " containing "
            << regionNCells[regioni] << " cells" << nl << endl;

        if (overwrite)
        {
            regionNames[regioni] = polyMesh::defaultRegion;
        }
        else if (regionCellZones[regioni] == -1)
        {
            regionNames[regioni] = "largest";
        }

        createAndWriteRegion
        (
            mesh,
            cellRegions,
            regionNames,
            faceInterfaces,
            interfacePatches,
            regioni,
            (overwrite ? oldInstance : runTime.name())
        );
    }
    else if (insidePoint)
    {
        const point insidePoint = args.optionRead<point>("insidePoint");

        const label celli = meshSearch::New(mesh).findCell(insidePoint);

        RemoteData<label> procCellRegion;
        if (celli != -1)
        {
            procCellRegion.proci = Pstream::myProcNo();
            procCellRegion.elementi = celli;
            procCellRegion.data = cellRegions[celli];
        }

        reduce(procCellRegion, RemoteData<label>::firstProcOp());

        if (procCellRegion.proci != -1)
        {
            Info<< "Found point " << insidePoint
                << " in cell " << procCellRegion.elementi;
            if (Pstream::parRun())
            {
                Info<< " of processor " << procCellRegion.proci;
            }
            Info<< nl << endl;
        }
        else
        {
            FatalErrorInFunction
                << "Did not find a cell containing the point " << insidePoint
                << ". The bounding box of the mesh is " << mesh.bounds()
                << "." << exit(FatalError);
        }

        const label regioni = procCellRegion.data;

        Info<< "Subsetting region " << regionNames[regioni]
            << " containing point " << insidePoint << nl << endl;

        if (overwrite)
        {
            regionNames[regioni] = polyMesh::defaultRegion;
        }
        else if (regionCellZones[regioni] == -1)
        {
            regionNames[regioni] = "insidePoint";
        }

        createAndWriteRegion
        (
            mesh,
            cellRegions,
            regionNames,
            faceInterfaces,
            interfacePatches,
            regioni,
            (overwrite ? oldInstance : runTime.name())
        );
    }
    else
    {
        forAll(regionNames, regioni)
        {
            createAndWriteRegion
            (
                mesh,
                cellRegions,
                regionNames,
                faceInterfaces,
                interfacePatches,
                regioni,
                (overwrite ? oldInstance : runTime.name())
            );
        }
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
