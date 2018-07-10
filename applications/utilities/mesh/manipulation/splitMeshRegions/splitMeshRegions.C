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
    splitMeshRegions

Description
    Splits mesh into multiple regions.

    Each region is defined as a domain whose cells can all be reached by
    cell-face-cell walking without crossing
    - boundary faces
    - additional faces from faceset (-blockedFaces faceSet).
    - any face in between differing cellZones (-cellZones)

    Output is:
    - volScalarField with regions as different scalars (-detectOnly)
            or
    - mesh with multiple regions and mapped patches. These patches
      either cover the whole interface between two region (default) or
      only part according to faceZones (-useFaceZones)
            or
    - mesh with cells put into cellZones (-makeCellZones)

    Note:
    - cellZonesOnly does not do a walk and uses the cellZones only. Use
    this if you don't mind having disconnected domains in a single region.
    This option requires all cells to be in one (and one only) cellZone.

    - cellZonesFileOnly behaves like -cellZonesOnly but reads the cellZones
    from the specified file. This allows one to explicitly specify the region
    distribution and still have multiple cellZones per region.

    - useCellZonesOnly does not do a walk and uses the cellZones only. Use
    this if you don't mind having disconnected domains in a single region.
    This option requires all cells to be in one (and one only) cellZone.

    - prefixRegion prefixes all normal patches with region name (interface
    (patches already have region name prefix)

    - Should work in parallel.
    cellZones can differ on either side of processor boundaries in which case
    the faces get moved from processor patch to directMapped patch. Not
    the faces get moved from processor patch to mapped patch. Not
    very well tested.

    - If a cell zone gets split into more than one region it can detect
    the largest matching region (-sloppyCellZones). This will accept any
    region that covers more than 50% of the zone. It has to be a subset
    so cannot have any cells in any other zone.

    - If explicitly a single region has been selected (-largestOnly or
      -insidePoint) its region name will be either
        - name of a cellZone it matches to or
        - "largestOnly" respectively "insidePoint" or
        - polyMesh::defaultRegion if additionally -overwrite
          (so it will overwrite the input mesh!)

    - writes maps like decomposePar back to original mesh:
        - pointRegionAddressing : for every point in this region the point in
        the original mesh
        - cellRegionAddressing  :   ,,      cell                ,,  cell    ,,
        - faceRegionAddressing  :   ,,      face                ,,  face in
        the original mesh + 'turning index'. For a face in the same orientation
        this is the original facelabel+1, for a turned face this is -facelabel-1
        - boundaryRegionAddressing : for every patch in this region the
        patch in the original mesh (or -1 if added patch)

\*---------------------------------------------------------------------------*/

#include "SortableList.H"
#include "argList.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "EdgeMap.H"
#include "syncTools.H"
#include "ReadFields.H"
#include "mappedWallPolyPatch.H"
#include "fvMeshTools.H"
#include "zeroGradientFvPatchFields.H"

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
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    forAll(patchesToRename, i)
    {
        label patchi = patchesToRename[i];
        polyPatch& pp = polyPatches[patchi];

        if (isA<coupledPolyPatch>(pp))
        {
            WarningInFunction
                << "Encountered coupled patch " << pp.name()
                << ". Will only rename the patch itself,"
                << " not any referred patches."
                << " This might have to be done by hand."
                << endl;
        }

        pp.name() = prefix + '_' + pp.name();
    }
    // Recalculate any demand driven data (e.g. group to name lookup)
    polyPatches.updateMesh();
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
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );
    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< "Mapping field " << fld.name() << endl;

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

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchi)
        {
            if (addedPatches.found(patchi))
            {
                tSubFld.ref().boundaryFieldRef()[patchi] ==
                    typename GeoField::value_type(Zero);
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
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );
    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< "Mapping field " << fld.name() << endl;

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

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchi)
        {
            if (addedPatches.found(patchi))
            {
                tSubFld.ref().boundaryFieldRef()[patchi] ==
                    typename GeoField::value_type(Zero);
            }
        }

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


// Get region-region interface name and sizes.
// Returns interfaces as straight list for looping in identical order.
void getInterfaceSizes
(
    const polyMesh& mesh,
    const bool useFaceZones,
    const labelList& cellRegion,
    const wordList& regionNames,

    edgeList& interfaces,
    List<Pair<word>>& interfaceNames,
    labelList& interfaceSizes,
    labelList& faceToInterface
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
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
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
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
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
    interfaceSizes.setSize(nInterfaces);
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
            interfaceSizes[nInterfaces] = infoIter();

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
    Pstream::scatter(interfaceSizes);
    Pstream::scatter(regionsToInterface);

    // Mark all inter-region faces.
    faceToInterface.setSize(mesh.nFaces(), -1);

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
}


// Create mesh for region.
autoPtr<mapPolyMesh> createRegionMesh
(
    const fvMesh& mesh,
    // Region info
    const labelList& cellRegion,
    const label regioni,
    const word& regionName,
    // Interface info
    const labelList& interfacePatches,
    const labelList& faceToInterface,

    autoPtr<fvMesh>& newMesh
)
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
        // if (!exists(io.objectPath()))
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
        // if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }


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
        label interfacei = faceToInterface[facei];

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

    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        newMesh,
        IOobject
        (
            regionName,
            mesh.time().timeName(),
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
    const bool prefixRegion,
    const labelList& faceToInterface,
    const labelList& interfacePatches,
    const label regioni,
    const word& newMeshInstance
)
{
    Info<< "Creating mesh for region " << regioni
        << ' ' << regionNames[regioni] << endl;

    autoPtr<fvMesh> newMesh;
    autoPtr<mapPolyMesh> map = createRegionMesh
    (
        mesh,
        cellRegion,
        regioni,
        regionNames[regioni],
        interfacePatches,
        faceToInterface,
        newMesh
    );


    // Make map of all added patches
    labelHashSet addedPatches(2*interfacePatches.size());
    forAll(interfacePatches, interfacei)
    {
        addedPatches.insert(interfacePatches[interfacei]);
        addedPatches.insert(interfacePatches[interfacei]+1);
    }


    Info<< "Mapping fields" << endl;

    // Map existing fields
    newMesh().updateMesh(map());

    // Add subsetted fields
    subsetVolFields<volScalarField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volVectorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volSymmTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );

    subsetSurfaceFields<surfaceScalarField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceVectorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceSymmTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );


    const polyBoundaryMesh& newPatches = newMesh().boundaryMesh();
    newPatches.checkParallelSync(true);

    // Delete empty patches
    // ~~~~~~~~~~~~~~~~~~~~

    // Create reordering list to move patches-to-be-deleted to end
    labelList oldToNew(newPatches.size(), -1);
    DynamicList<label> sharedPatches(newPatches.size());
    label newI = 0;

    Info<< "Deleting empty patches" << endl;

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

    // Move all deleteable patches to the end
    forAll(oldToNew, patchi)
    {
        if (oldToNew[patchi] == -1)
        {
            oldToNew[patchi] = newI++;
        }
    }

    // reorderPatches(newMesh(), oldToNew, nNewPatches);
    fvMeshTools::reorderPatches(newMesh(), oldToNew, nNewPatches, true);

    // Rename shared patches with region name
    if (prefixRegion)
    {
        Info<< "Prefixing patches with region name" << endl;

        renamePatches(newMesh(), regionNames[regioni], sharedPatches);
    }


    Info<< "Writing new mesh" << endl;

    newMesh().setInstance(newMeshInstance);
    newMesh().write();

    // Write addressing files like decomposePar
    Info<< "Writing addressing to base mesh" << endl;

    labelIOList pointProcAddressing
    (
        IOobject
        (
            "pointRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        map().pointMap()
    );
    Info<< "Writing map " << pointProcAddressing.name()
        << " from region" << regioni
        << " points back to base mesh." << endl;
    pointProcAddressing.write();

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newMesh().nFaces()
    );
    forAll(faceProcAddressing, facei)
    {
        // face + turning index. (see decomposePar)
        // Is the face pointing in the same direction?
        label oldFacei = map().faceMap()[facei];

        if
        (
            map().cellMap()[newMesh().faceOwner()[facei]]
         == mesh.faceOwner()[oldFacei]
        )
        {
            faceProcAddressing[facei] = oldFacei+1;
        }
        else
        {
            faceProcAddressing[facei] = -(oldFacei+1);
        }
    }
    Info<< "Writing map " << faceProcAddressing.name()
        << " from region" << regioni
        << " faces back to base mesh." << endl;
    faceProcAddressing.write();

    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        map().cellMap()
    );
    Info<< "Writing map " <<cellProcAddressing.name()
        << " from region" << regioni
        << " cells back to base mesh." << endl;
    cellProcAddressing.write();

    labelIOList boundaryProcAddressing
    (
        IOobject
        (
            "boundaryRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(nNewPatches, -1)
    );
    forAll(oldToNew, i)
    {
        if (!addedPatches.found(i))
        {
            label newI = oldToNew[i];
            if (newI >= 0 && newI < nNewPatches)
            {
                boundaryProcAddressing[oldToNew[i]] = i;
            }
        }
    }
    Info<< "Writing map " << boundaryProcAddressing.name()
        << " from region" << regioni
        << " boundary back to base mesh." << endl;
    boundaryProcAddressing.write();
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
    Info<< nl << "Adding patches" << nl << endl;

    labelList interfacePatches(interfaces.size());

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];
        const Pair<word>& names = interfaceNames[interI];

        // Info<< "For interface " << interI
        //    << " between regions " << e
        //    << " trying to add patches " << names << endl;


        mappedWallPolyPatch patch1
        (
            names[0],
            0,                  // overridden
            0,                  // overridden
            0,                  // overridden
            regionNames[e[1]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[1],           // samplePatch
            point::zero,        // offset
            mesh.boundaryMesh()
        );

        interfacePatches[interI] = fvMeshTools::addPatch
        (
            mesh,
            patch1,
            dictionary(),   // optional per field value
            calculatedFvPatchField<scalar>::typeName,
            true            // validBoundary
        );

        mappedWallPolyPatch patch2
        (
            names[1],
            0,
            0,
            0,
            regionNames[e[0]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[0],
            point::zero,        // offset
            mesh.boundaryMesh()
        );
        fvMeshTools::addPatch
        (
            mesh,
            patch2,
            dictionary(),   // optional per field value
            calculatedFvPatchField<scalar>::typeName,
            true            // validBoundary
        );

        Info<< "For interface between region " << regionNames[e[0]]
            << " and " << regionNames[e[1]] << " added patches" << endl
            << "    " << interfacePatches[interI]
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]].name()
            << endl
            << "    " << interfacePatches[interI]+1
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]+1].name()
            << endl;
    }
    return interfacePatches;
}


// Find region that covers most of cell zone
label findCorrespondingRegion
(
    const labelList& existingZoneID,    // per cell the (unique) zoneID
    const labelList& cellRegion,
    const label nCellRegions,
    const label zoneI,
    const label minOverlapSize
)
{
    // Per region the number of cells in zoneI
    labelList cellsInZone(nCellRegions, 0);

    forAll(cellRegion, celli)
    {
        if (existingZoneID[celli] == zoneI)
        {
            cellsInZone[cellRegion[celli]]++;
        }
    }

    Pstream::listCombineGather(cellsInZone, plusEqOp<label>());
    Pstream::listCombineScatter(cellsInZone);

    // Pick region with largest overlap of zoneI
    label regioni = findMax(cellsInZone);


    if (cellsInZone[regioni] < minOverlapSize)
    {
        // Region covers too little of zone. Not good enough.
        regioni = -1;
    }
    else
    {
        // Check that region contains no cells that aren't in cellZone.
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == regioni && existingZoneID[celli] != zoneI)
            {
                // celli in regioni but not in zoneI
                regioni = -1;
                break;
            }
        }
        // If one in error, all should be in error. Note that branch gets taken
        // on all procs.
        reduce(regioni, minOp<label>());
    }

    return regioni;
}


// Get zone per cell
// - non-unique zoning
// - coupled zones
void getZoneID
(
    const polyMesh& mesh,
    const cellZoneMesh& cellZones,
    labelList& zoneID,
    labelList& neiZoneID
)
{
    // Existing zoneID
    zoneID.setSize(mesh.nCells());
    zoneID = -1;

    forAll(cellZones, zoneI)
    {
        const cellZone& cz = cellZones[zoneI];

        forAll(cz, i)
        {
            label celli = cz[i];
            if (zoneID[celli] == -1)
            {
                zoneID[celli] = zoneI;
            }
            else
            {
                FatalErrorInFunction
                    << "Cell " << celli << " with cell centre "
                    << mesh.cellCentres()[celli]
                    << " is multiple zones. This is not allowed." << endl
                    << "It is in zone " << cellZones[zoneID[celli]].name()
                    << " and in zone " << cellZones[zoneI].name()
                    << exit(FatalError);
            }
        }
    }

    // Neighbour zoneID.
    neiZoneID.setSize(mesh.nFaces()-mesh.nInternalFaces());

    forAll(neiZoneID, i)
    {
        neiZoneID[i] = zoneID[mesh.faceOwner()[i+mesh.nInternalFaces()]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiZoneID);
}


void matchRegions
(
    const bool sloppyCellZones,
    const polyMesh& mesh,

    const label nCellRegions,
    const labelList& cellRegion,

    const word& defaultRegionName,

    labelList& regionToZone,
    wordList& regionNames,
    labelList& zoneToRegion
)
{
    const cellZoneMesh& cellZones = mesh.cellZones();

    regionToZone.setSize(nCellRegions, -1);
    regionNames.setSize(nCellRegions);
    zoneToRegion.setSize(cellZones.size(), -1);

    // Get current per cell zoneID
    labelList zoneID(mesh.nCells(), -1);
    labelList neiZoneID(mesh.nFaces()-mesh.nInternalFaces());
    getZoneID(mesh, cellZones, zoneID, neiZoneID);

    // Sizes per cellzone
    labelList zoneSizes(cellZones.size(), 0);
    {
        List<wordList> zoneNames(Pstream::nProcs());
        zoneNames[Pstream::myProcNo()] = cellZones.names();
        Pstream::gatherList(zoneNames);
        Pstream::scatterList(zoneNames);

        forAll(zoneNames, proci)
        {
            if (zoneNames[proci] != zoneNames[0])
            {
                FatalErrorInFunction
                    << "cellZones not synchronised across processors." << endl
                    << "Master has cellZones " << zoneNames[0] << endl
                    << "Processor " << proci
                    << " has cellZones " << zoneNames[proci]
                    << exit(FatalError);
            }
        }

        forAll(cellZones, zoneI)
        {
            zoneSizes[zoneI] = returnReduce
            (
                cellZones[zoneI].size(),
                sumOp<label>()
            );
        }
    }


    if (sloppyCellZones)
    {
        Info<< "Trying to match regions to existing cell zones;"
            << " region can be subset of cell zone." << nl << endl;

        forAll(cellZones, zoneI)
        {
            label regioni = findCorrespondingRegion
            (
                zoneID,
                cellRegion,
                nCellRegions,
                zoneI,
                label(0.5*zoneSizes[zoneI]) // minimum overlap
            );

            if (regioni != -1)
            {
                Info<< "Sloppily matched region " << regioni
                    //<< " size " << regionSizes[regioni]
                    << " to zone " << zoneI << " size " << zoneSizes[zoneI]
                    << endl;
                zoneToRegion[zoneI] = regioni;
                regionToZone[regioni] = zoneI;
                regionNames[regioni] = cellZones[zoneI].name();
            }
        }
    }
    else
    {
        Info<< "Trying to match regions to existing cell zones." << nl << endl;

        forAll(cellZones, zoneI)
        {
            label regioni = findCorrespondingRegion
            (
                zoneID,
                cellRegion,
                nCellRegions,
                zoneI,
                1               // minimum overlap
            );

            if (regioni != -1)
            {
                zoneToRegion[zoneI] = regioni;
                regionToZone[regioni] = zoneI;
                regionNames[regioni] = cellZones[zoneI].name();
            }
        }
    }

    label nUnmatchedRegions = 0;

    forAll(regionToZone, regioni)
    {
        if (regionToZone[regioni] == -1)
        {
            nUnmatchedRegions++;
        }
    }

    if (nUnmatchedRegions)
    {
        label nUnmatchedi = 1;

        // Allocate region names for unmatched regions.
        forAll(regionToZone, regioni)
        {
            if (regionToZone[regioni] == -1)
            {
                if
                (
                    nUnmatchedRegions == 1
                 && defaultRegionName != standardRegionName
                )
                {
                    regionNames[regioni] = defaultRegionName;
                }
                else
                {
                    regionNames[regioni] =
                        defaultRegionName + Foam::name(nUnmatchedi++);
                }
            }
        }
    }
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

        Info<< "Writing region per cell file (for manual decomposition) to "
            << cellToRegion.objectPath() << nl << endl;
    }
    // Write for postprocessing
    {
        volScalarField cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(cellRegion, celli)
        {
            cellToRegion[celli] = cellRegion[celli];
        }
        cellToRegion.write();

        Info<< "Writing region per cell as volScalarField to "
            << cellToRegion.objectPath() << nl << endl;
    }
}




int main(int argc, char *argv[])
{
    argList::addNote
    (
        "splits mesh into multiple regions (detected by walking across faces)"
    );
    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    argList::addBoolOption
    (
        "cellZones",
        "additionally split cellZones off into separate regions"
    );
    argList::addBoolOption
    (
        "cellZonesOnly",
        "use cellZones only to split mesh into regions; do not use walking"
    );
    argList::addOption
    (
        "cellZonesFileOnly",
        "file",
        "like -cellZonesOnly, but use specified file"
    );
    argList::addOption
    (
        "blockedFaces",
        "faceSet",
        "specify additional region boundaries that walking does not cross"
    );
    argList::addBoolOption
    (
        "makeCellZones",
        "place cells into cellZones instead of splitting mesh"
    );
    argList::addBoolOption
    (
        "largestOnly",
        "only write largest region"
    );
    argList::addOption
    (
        "insidePoint",
        "point",
        "only write region containing point"
    );
    argList::addBoolOption
    (
        "detectOnly",
        "do not write mesh"
    );
    argList::addBoolOption
    (
        "sloppyCellZones",
        "try to match heuristically regions to existing cell zones"
    );
    argList::addBoolOption
    (
        "useFaceZones",
        "use faceZones to patch inter-region faces instead of single patch"
    );
    argList::addBoolOption
    (
        "prefixRegion",
        "prefix region name to all patches, not just coupling patches"
    );
    argList::addOption
    (
        "defaultRegionName",
        "name",
        "base name of the unspecified regions, defaults to \"region\""
    );

    #include "setRootCase.H"
    #include "createTime.H"

    runTime.functionObjects().off();

    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    word blockedFacesName;
    if (args.optionReadIfPresent("blockedFaces", blockedFacesName))
    {
        Info<< "Reading blocked internal faces from faceSet "
            << blockedFacesName << nl << endl;
    }

    const bool makeCellZones    = args.optionFound("makeCellZones");
    const bool largestOnly      = args.optionFound("largestOnly");
    const bool insidePoint      = args.optionFound("insidePoint");
    const bool useCellZones     = args.optionFound("cellZones");
    const bool useCellZonesOnly = args.optionFound("cellZonesOnly");
    const bool useCellZonesFile = args.optionFound("cellZonesFileOnly");
    const bool overwrite        = args.optionFound("overwrite");
    const bool detectOnly       = args.optionFound("detectOnly");
    const bool sloppyCellZones  = args.optionFound("sloppyCellZones");
    const bool useFaceZones     = args.optionFound("useFaceZones");
    const bool prefixRegion     = args.optionFound("prefixRegion");

    if
    (
        (useCellZonesOnly || useCellZonesFile)
     && (useCellZones || blockedFacesName.size())
    )
    {
        FatalErrorInFunction
            << "You cannot specify both -cellZonesOnly or -cellZonesFileOnly"
            << " (which specify complete split)"
            << " in combination with -blockedFaces or -cellZones"
            << " (which imply a split based on topology)"
            << exit(FatalError);
    }


    if (useFaceZones)
    {
        Info<< "Using current faceZones to divide inter-region interfaces"
            << " into multiple patches."
            << nl << endl;
    }
    else
    {
        Info<< "Creating single patch per inter-region interface."
            << nl << endl;
    }


    if (insidePoint && largestOnly)
    {
        FatalErrorInFunction
            << "You cannot specify both -largestOnly"
            << " (keep region with most cells)"
            << " and -insidePoint (keep region containing point)"
            << exit(FatalError);
    }

    const word defaultRegionName
    (
        args.optionLookupOrDefault("defaultRegionName", standardRegionName)
    );

    const cellZoneMesh& cellZones = mesh.cellZones();

    // Existing zoneID
    labelList zoneID(mesh.nCells(), -1);
    // Neighbour zoneID.
    labelList neiZoneID(mesh.nFaces()-mesh.nInternalFaces());
    getZoneID(mesh, cellZones, zoneID, neiZoneID);



    // Determine per cell the region it belongs to
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // cellRegion is the labelList with the region per cell.
    labelList cellRegion;
    // Region per zone
    labelList regionToZone;
    // Name of region
    wordList regionNames;
    // Zone to region
    labelList zoneToRegion;

    label nCellRegions = 0;
    if (useCellZonesOnly)
    {
        Info<< "Using current cellZones to split mesh into regions."
            << " This requires all"
            << " cells to be in one and only one cellZone." << nl << endl;

        label unzonedCelli = findIndex(zoneID, -1);
        if (unzonedCelli != -1)
        {
            FatalErrorInFunction
                << "For the cellZonesOnly option all cells "
                << "have to be in a cellZone." << endl
                << "Cell " << unzonedCelli
                << " at" << mesh.cellCentres()[unzonedCelli]
                << " is not in a cellZone. There might be more unzoned cells."
                << exit(FatalError);
        }
        cellRegion = zoneID;
        nCellRegions = gMax(cellRegion)+1;
        regionToZone.setSize(nCellRegions);
        regionNames.setSize(nCellRegions);
        zoneToRegion.setSize(cellZones.size(), -1);
        for (label regioni = 0; regioni < nCellRegions; regioni++)
        {
            regionToZone[regioni] = regioni;
            zoneToRegion[regioni] = regioni;
            regionNames[regioni] = cellZones[regioni].name();
        }
    }
    else if (useCellZonesFile)
    {
        const word zoneFile = args.option("cellZonesFileOnly");
        Info<< "Reading split from cellZones file " << zoneFile << endl
            << "This requires all"
            << " cells to be in one and only one cellZone." << nl << endl;

        cellZoneMesh newCellZones
        (
            IOobject
            (
                zoneFile,
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh
        );

        labelList newZoneID(mesh.nCells(), -1);
        labelList newNeiZoneID(mesh.nFaces()-mesh.nInternalFaces());
        getZoneID(mesh, newCellZones, newZoneID, newNeiZoneID);

        label unzonedCelli = findIndex(newZoneID, -1);
        if (unzonedCelli != -1)
        {
            FatalErrorInFunction
                << "For the cellZonesFileOnly option all cells "
                << "have to be in a cellZone." << endl
                << "Cell " << unzonedCelli
                << " at" << mesh.cellCentres()[unzonedCelli]
                << " is not in a cellZone. There might be more unzoned cells."
                << exit(FatalError);
        }
        cellRegion = newZoneID;
        nCellRegions = gMax(cellRegion)+1;
        zoneToRegion.setSize(newCellZones.size(), -1);
        regionToZone.setSize(nCellRegions);
        regionNames.setSize(nCellRegions);
        for (label regioni = 0; regioni < nCellRegions; regioni++)
        {
            regionToZone[regioni] = regioni;
            zoneToRegion[regioni] = regioni;
            regionNames[regioni] = newCellZones[regioni].name();
        }
    }
    else
    {
        // Determine connected regions
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Mark additional faces that are blocked
        boolList blockedFace;

        // Read from faceSet
        if (blockedFacesName.size())
        {
            faceSet blockedFaceSet(mesh, blockedFacesName);
            Info<< "Read "
                << returnReduce(blockedFaceSet.size(), sumOp<label>())
                << " blocked faces from set " << blockedFacesName << nl << endl;

            blockedFace.setSize(mesh.nFaces(), false);

            forAllConstIter(faceSet, blockedFaceSet, iter)
            {
                blockedFace[iter.key()] = true;
            }
        }

        // Imply from differing cellZones
        if (useCellZones)
        {
            blockedFace.setSize(mesh.nFaces(), false);

            for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
            {
                label own = mesh.faceOwner()[facei];
                label nei = mesh.faceNeighbour()[facei];

                if (zoneID[own] != zoneID[nei])
                {
                    blockedFace[facei] = true;
                }
            }

            // Different cellZones on either side of processor patch.
            forAll(neiZoneID, i)
            {
                label facei = i+mesh.nInternalFaces();

                if (zoneID[mesh.faceOwner()[facei]] != neiZoneID[i])
                {
                    blockedFace[facei] = true;
                }
            }
        }

        // Do a topological walk to determine regions
        regionSplit regions(mesh, blockedFace);
        nCellRegions = regions.nRegions();
        cellRegion.transfer(regions);

        // Make up region names. If possible match them to existing zones.
        matchRegions
        (
            sloppyCellZones,
            mesh,
            nCellRegions,
            cellRegion,
            defaultRegionName,

            regionToZone,
            regionNames,
            zoneToRegion
        );

        // Override any default region names if single region selected
        if (largestOnly || insidePoint)
        {
            forAll(regionToZone, regioni)
            {
                if (regionToZone[regioni] == -1)
                {
                    if (overwrite)
                    {
                        regionNames[regioni] = polyMesh::defaultRegion;
                    }
                    else if (insidePoint)
                    {
                        regionNames[regioni] = "insidePoint";
                    }
                    else if (largestOnly)
                    {
                        regionNames[regioni] = "largestOnly";
                    }
                }
            }
        }
    }

    Info<< endl << "Number of regions:" << nCellRegions << nl << endl;


    // Write decomposition to file
    writeCellToRegion(mesh, cellRegion);



    // Sizes per region
    // ~~~~~~~~~~~~~~~~

    labelList regionSizes(nCellRegions, 0);

    forAll(cellRegion, celli)
    {
        regionSizes[cellRegion[celli]]++;
    }
    forAll(regionSizes, regioni)
    {
        reduce(regionSizes[regioni], sumOp<label>());
    }

    Info<< "Region\tCells" << nl
        << "------\t-----" << endl;

    forAll(regionSizes, regioni)
    {
        Info<< regioni << '\t' << regionSizes[regioni] << nl;
    }
    Info<< endl;



    // Print region to zone
    Info<< "Region\tZone\tName" << nl
        << "------\t----\t----" << endl;
    forAll(regionToZone, regioni)
    {
        Info<< regioni << '\t' << regionToZone[regioni] << '\t'
            << regionNames[regioni] << nl;
    }
    Info<< endl;



    // Since we're going to mess with patches and zones make sure all
    // is synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    mesh.faceZones().checkParallelSync(true);


    // Interfaces
    // ----------
    // per interface:
    // - the two regions on either side
    // - the name
    // - the (global) size
    edgeList interfaces;
    List<Pair<word>> interfaceNames;
    labelList interfaceSizes;
    // per face the interface
    labelList faceToInterface;

    getInterfaceSizes
    (
        mesh,
        useFaceZones,
        cellRegion,
        regionNames,

        interfaces,
        interfaceNames,
        interfaceSizes,
        faceToInterface
    );

    Info<< "Sizes of interfaces between regions:" << nl << nl
        << "Interface\tRegion\tRegion\tFaces" << nl
        << "---------\t------\t------\t-----" << endl;

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];

        Info<< interI
            << "\t\t" << e[0] << '\t' << e[1]
            << '\t' << interfaceSizes[interI] << nl;
    }
    Info<< endl;


    if (detectOnly)
    {
        return 0;
    }


    // Read objects in time directory
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);

    Info<< endl;


    // Remove any demand-driven fields ('S', 'V' etc)
    mesh.clearOut();


    if (nCellRegions == 1)
    {
        Info<< "Only one region. Doing nothing." << endl;
    }
    else if (makeCellZones)
    {
        Info<< "Putting cells into cellZones instead of splitting mesh."
            << endl;

        // Check if region overlaps with existing zone. If so keep.

        for (label regioni = 0; regioni < nCellRegions; regioni++)
        {
            label zoneI = regionToZone[regioni];

            if (zoneI != -1)
            {
                Info<< "    Region " << regioni << " : corresponds to existing"
                    << " cellZone "
                    << zoneI << ' ' << cellZones[zoneI].name() << endl;
            }
            else
            {
                // Create new cellZone.
                labelList regionCells = findIndices(cellRegion, regioni);

                word zoneName = "region" + Foam::name(regioni);

                zoneI = cellZones.findZoneID(zoneName);

                if (zoneI == -1)
                {
                    zoneI = cellZones.size();
                    mesh.cellZones().setSize(zoneI+1);
                    mesh.cellZones().set
                    (
                        zoneI,
                        new cellZone
                        (
                            zoneName,           // name
                            regionCells,        // addressing
                            zoneI,              // index
                            cellZones           // cellZoneMesh
                        )
                    );
                }
                else
                {
                    mesh.cellZones()[zoneI].clearAddressing();
                    mesh.cellZones()[zoneI] = regionCells;
                }
                Info<< "    Region " << regioni << " : created new cellZone "
                    << zoneI << ' ' << cellZones[zoneI].name() << endl;
            }
        }
        mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;

        if (!overwrite)
        {
            runTime++;
            mesh.setInstance(runTime.timeName());
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        Info<< "Writing cellZones as new mesh to time " << runTime.timeName()
            << nl << endl;

        mesh.write();


        // Write cellSets for convenience
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Info<< "Writing cellSets corresponding to cellZones." << nl << endl;

        forAll(cellZones, zoneI)
        {
            const cellZone& cz = cellZones[zoneI];

            cellSet(mesh, cz.name(), cz).write();
        }
    }
    else
    {
        // Add patches for interfaces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Add all possible patches. Empty ones get filtered later on.
        Info<< nl << "Adding patches" << nl << endl;

        labelList interfacePatches
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
        // ~~~~~~~~~~~~~~

        if (insidePoint)
        {
            const point insidePoint = args.optionRead<point>("insidePoint");

            label regioni = -1;

            (void)mesh.tetBasePtIs();

            label celli = mesh.findCell(insidePoint);

            Info<< nl << "Found point " << insidePoint << " in cell " << celli
                << endl;

            if (celli != -1)
            {
                regioni = cellRegion[celli];
            }

            reduce(regioni, maxOp<label>());

            Info<< nl
                << "Subsetting region " << regioni
                << " containing point " << insidePoint << endl;

            if (regioni == -1)
            {
                FatalErrorInFunction
                    << "Point " << insidePoint
                    << " is not inside the mesh." << nl
                    << "Bounding box of the mesh:" << mesh.bounds()
                    << exit(FatalError);
            }

            createAndWriteRegion
            (
                mesh,
                cellRegion,
                regionNames,
                prefixRegion,
                faceToInterface,
                interfacePatches,
                regioni,
                (overwrite ? oldInstance : runTime.timeName())
            );
        }
        else if (largestOnly)
        {
            label regioni = findMax(regionSizes);

            Info<< nl
                << "Subsetting region " << regioni
                << " of size " << regionSizes[regioni]
                << " as named region " << regionNames[regioni] << endl;

            createAndWriteRegion
            (
                mesh,
                cellRegion,
                regionNames,
                prefixRegion,
                faceToInterface,
                interfacePatches,
                regioni,
                (overwrite ? oldInstance : runTime.timeName())
            );
        }
        else
        {
            // Split all
            for (label regioni = 0; regioni < nCellRegions; regioni++)
            {
                Info<< nl
                    << "Region " << regioni << nl
                    << "-------- " << endl;

                createAndWriteRegion
                (
                    mesh,
                    cellRegion,
                    regionNames,
                    prefixRegion,
                    faceToInterface,
                    interfacePatches,
                    regioni,
                    (overwrite ? oldInstance : runTime.timeName())
                );
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
