/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    renumberMesh

Description
    Renumbers the cell list in order to reduce the bandwidth, reading and
    renumbering all fields from all the time directories.

    By default uses bandCompression (CuthillMcKee) but will
    read system/renumberMeshDict if -dict option is present

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "polyTopoChange.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SortableList.H"
#include "decompositionMethod.H"
#include "renumberMethod.H"
#include "zeroGradientFvPatchFields.H"
#include "CuthillMcKeeRenumber.H"
#include "fvMeshSubset.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "systemDict.H"

#ifdef FOAM_USE_ZOLTAN
    #include "zoltanRenumber.H"
#endif


using namespace Foam;


// Create named field from labelList for postprocessing
tmp<volScalarField> createScalarField
(
    const fvMesh& mesh,
    const word& name,
    const labelList& elems
)
{
    tmp<volScalarField> tfld
    (
        volScalarField::New
        (
            name,
            mesh,
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& fld = tfld.ref();

    forAll(fld, celli)
    {
       fld[celli] = elems[celli];
    }

    return tfld;
}


// Calculate band of matrix
label getBand(const labelList& owner, const labelList& neighbour)
{
    label band = 0;

    forAll(neighbour, facei)
    {
        label diff = neighbour[facei] - owner[facei];

        if (diff > band)
        {
            band = diff;
        }
    }
    return band;
}


// Calculate band of matrix
void getBand
(
    const bool calculateIntersect,
    const label nCells,
    const labelList& owner,
    const labelList& neighbour,
    label& bandwidth,
    scalar& profile,            // scalar to avoid overflow
    scalar& sumSqrIntersect     // scalar to avoid overflow
)
{
    labelList cellBandwidth(nCells, 0);
    scalarField nIntersect(nCells, 0.0);

    forAll(neighbour, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // Note: mag not necessary for correct (upper-triangular) ordering.
        label diff = nei-own;
        cellBandwidth[nei] = max(cellBandwidth[nei], diff);
    }

    bandwidth = max(cellBandwidth);

    // Do not use field algebra because of conversion label to scalar
    profile = 0.0;
    forAll(cellBandwidth, celli)
    {
        profile += 1.0*cellBandwidth[celli];
    }

    sumSqrIntersect = 0.0;
    if (calculateIntersect)
    {
        forAll(nIntersect, celli)
        {
            for (label colI = celli-cellBandwidth[celli]; colI <= celli; colI++)
            {
                nIntersect[colI] += 1.0;
            }
        }

        sumSqrIntersect = sum(Foam::sqr(nIntersect));
    }
}


// Determine upper-triangular face order
labelList getFaceOrder
(
    const primitiveMesh& mesh,
    const labelList& cellOrder      // New to old cell
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));

    labelList oldToNewFace(mesh.nFaces(), -1);

    label newFacei = 0;

    labelList nbr;
    labelList order;

    forAll(cellOrder, newCelli)
    {
        label oldCelli = cellOrder[newCelli];

        const cell& cFaces = mesh.cells()[oldCelli];

        // Neighbouring cells
        nbr.setSize(cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh.isInternalFace(facei))
            {
                // Internal face. Get cell on other side.
                label nbrCelli = reverseCellOrder[mesh.faceNeighbour()[facei]];
                if (nbrCelli == newCelli)
                {
                    nbrCelli = reverseCellOrder[mesh.faceOwner()[facei]];
                }

                if (newCelli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        order.setSize(nbr.size());
        sortedOrder(nbr, order);

        forAll(order, i)
        {
            label index = order[i];
            if (nbr[index] != -1)
            {
                oldToNewFace[cFaces[index]] = newFacei++;
            }
        }
    }

    // Leave patch faces intact.
    for (label facei = newFacei; facei < mesh.nFaces(); facei++)
    {
        oldToNewFace[facei] = facei;
    }


    // Check done all faces.
    forAll(oldToNewFace, facei)
    {
        if (oldToNewFace[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position" << " for face " << facei
                << abort(FatalError);
        }
    }

    return invert(mesh.nFaces(), oldToNewFace);
}


// Determine face order such that inside region faces are sorted
// upper-triangular but in between region faces are handled like boundary faces.
labelList getRegionFaceOrder
(
    const primitiveMesh& mesh,
    const labelList& cellOrder,     // New to old cell
    const labelList& cellToRegion   // Old cell to region
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));

    labelList oldToNewFace(mesh.nFaces(), -1);

    label newFacei = 0;

    label prevRegion = -1;

    forAll(cellOrder, newCelli)
    {
        label oldCelli = cellOrder[newCelli];

        if (cellToRegion[oldCelli] != prevRegion)
        {
            prevRegion = cellToRegion[oldCelli];
        }

        const cell& cFaces = mesh.cells()[oldCelli];

        SortableList<label> nbr(cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh.isInternalFace(facei))
            {
                // Internal face. Get cell on other side.
                label nbrCelli = reverseCellOrder[mesh.faceNeighbour()[facei]];
                if (nbrCelli == newCelli)
                {
                    nbrCelli = reverseCellOrder[mesh.faceOwner()[facei]];
                }

                if (cellToRegion[oldCelli] != cellToRegion[cellOrder[nbrCelli]])
                {
                    // Treat like external face. Do later.
                    nbr[i] = -1;
                }
                else if (newCelli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNewFace[cFaces[nbr.indices()[i]]] = newFacei++;
            }
        }
    }

    // Do region interfaces
    label nRegions = max(cellToRegion)+1;
    {
        // Sort in increasing region
        SortableList<label> sortKey(mesh.nFaces(), labelMax);

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label ownRegion = cellToRegion[mesh.faceOwner()[facei]];
            label neiRegion = cellToRegion[mesh.faceNeighbour()[facei]];

            if (ownRegion != neiRegion)
            {
                sortKey[facei] =
                    min(ownRegion, neiRegion)*nRegions
                   +max(ownRegion, neiRegion);
            }
        }
        sortKey.sort();

        // Extract.
        label prevKey = -1;
        forAll(sortKey, i)
        {
            label key = sortKey[i];

            if (key == labelMax)
            {
                break;
            }

            if (prevKey != key)
            {
                prevKey = key;
            }

            oldToNewFace[sortKey.indices()[i]] = newFacei++;
        }
    }

    // Leave patch faces intact.
    for (label facei = newFacei; facei < mesh.nFaces(); facei++)
    {
        oldToNewFace[facei] = facei;
    }


    // Check done all faces.
    forAll(oldToNewFace, facei)
    {
        if (oldToNewFace[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position"
                << " for face " << facei
                << abort(FatalError);
        }
    }

    return invert(mesh.nFaces(), oldToNewFace);
}


// cellOrder: old cell for every new cell
// faceOrder: old face for every new face. Ordering of boundary faces not
//     changed.
autoPtr<mapPolyMesh> reorderMesh
(
    polyMesh& mesh,
    const labelList& cellOrder,
    const labelList& faceOrder
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));
    labelList reverseFaceOrder(invert(faceOrder.size(), faceOrder));

    faceList newFaces(reorder(reverseFaceOrder, mesh.faces()));
    labelList newOwner
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceOwner())
        )
    );
    labelList newNeighbour
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceNeighbour())
        )
    );

    // Check if any faces need swapping.
    labelHashSet flipFaceFlux(newOwner.size());
    forAll(newNeighbour, facei)
    {
        label own = newOwner[facei];
        label nei = newNeighbour[facei];

        if (nei < own)
        {
            newFaces[facei].flip();
            Swap(newOwner[facei], newNeighbour[facei]);
            flipFaceFlux.insert(facei);
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelList patchSizes(patches.size());
    labelList patchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    labelListList patchPointMap(patches.size());

    forAll(patches, patchi)
    {
        patchSizes[patchi] = patches[patchi].size();
        patchStarts[patchi] = patches[patchi].start();
        oldPatchNMeshPoints[patchi] = patches[patchi].nPoints();
        patchPointMap[patchi] = identity(patches[patchi].nPoints());
    }

    mesh.resetPrimitives
    (
        NullObjectMove<pointField>(),
        move(newFaces),
        move(newOwner),
        move(newNeighbour),
        patchSizes,
        patchStarts,
        true
    );


    // Re-do the faceZones
    {
        meshFaceZones& faceZones = mesh.faceZones();
        faceZones.clearAddressing();
        forAll(faceZones, zoneI)
        {
            faceZone& fZone = faceZones[zoneI];
            labelList newAddressing(fZone.size());
            boolList newFlipMap(fZone.size());
            forAll(fZone, i)
            {
                label oldFacei = fZone[i];
                newAddressing[i] = reverseFaceOrder[oldFacei];
                if (flipFaceFlux.found(newAddressing[i]))
                {
                    newFlipMap[i] = !fZone.flipMap()[i];
                }
                else
                {
                    newFlipMap[i] = fZone.flipMap()[i];
                }
            }
            labelList newToOld;
            sortedOrder(newAddressing, newToOld);
            fZone.resetAddressing
            (
                UIndirectList<label>(newAddressing, newToOld)(),
                UIndirectList<bool>(newFlipMap, newToOld)()
            );
        }
    }
    // Re-do the cellZones
    {
        meshCellZones& cellZones = mesh.cellZones();
        cellZones.clearAddressing();
        forAll(cellZones, zoneI)
        {
            cellZones[zoneI] = UIndirectList<label>
            (
                reverseCellOrder,
                cellZones[zoneI]
            )();
            Foam::sort(cellZones[zoneI]);
        }
    }


    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            mesh,                       // const polyMesh& mesh,
            mesh.nPoints(),             // nOldPoints,
            mesh.nFaces(),              // nOldFaces,
            mesh.nCells(),              // nOldCells,
            identity(mesh.nPoints()),   // pointMap,
            List<objectMap>(0),         // pointsFromPoints,
            faceOrder,                  // faceMap,
            List<objectMap>(0),         // facesFromPoints,
            List<objectMap>(0),         // facesFromEdges,
            List<objectMap>(0),         // facesFromFaces,
            cellOrder,                  // cellMap,
            List<objectMap>(0),         // cellsFromPoints,
            List<objectMap>(0),         // cellsFromEdges,
            List<objectMap>(0),         // cellsFromFaces,
            List<objectMap>(0),         // cellsFromCells,
            identity(mesh.nPoints()),   // reversePointMap,
            reverseFaceOrder,           // reverseFaceMap,
            reverseCellOrder,           // reverseCellMap,
            flipFaceFlux,               // flipFaceFlux,
            patchPointMap,              // patchPointMap,
            labelListList(0),           // pointZoneMap,
            labelListList(0),           // faceZonePointMap,
            labelListList(0),           // faceZoneFaceMap,
            labelListList(0),           // cellZoneMap,
            pointField(0),              // preMotionPoints,
            patchStarts,                // oldPatchStarts,
            oldPatchNMeshPoints,        // oldPatchNMeshPoints
            autoPtr<scalarField>()      // oldCellVolumes
        )
    );
}


// Return new to old cell numbering
labelList regionRenumber
(
    const renumberMethod& method,
    const fvMesh& mesh,
    const labelList& cellToRegion
)
{
    Info<< "Determining cell order:" << endl;

    labelList cellOrder(cellToRegion.size());

    label nRegions = max(cellToRegion)+1;

    labelListList regionToCells(invertOneToMany(nRegions, cellToRegion));

    label celli = 0;

    forAll(regionToCells, regionI)
    {
        Info<< "    region " << regionI << " starts at " << celli << endl;

        // Make sure no parallel comms
        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        // Per region do a reordering.
        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(cellToRegion, regionI);

        const fvMesh& subMesh = subsetter.subMesh();

        labelList subCellOrder = method.renumber
        (
            subMesh,
            subMesh.cellCentres()
        );

        // Restore state
        UPstream::parRun() = oldParRun;

        const labelList& cellMap = subsetter.cellMap();

        forAll(subCellOrder, i)
        {
            cellOrder[celli++] = cellMap[subCellOrder[i]];
        }
    }
    Info<< endl;

    return cellOrder;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Renumber mesh to minimise bandwidth"
    );

    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    #include "addTimeOptions.H"
    #include "addDictOption.H"
    argList::addBoolOption
    (
        "frontWidth",
        "calculate the rms of the frontwidth"
    );
    argList::addBoolOption
    (
        "noFields",
        "do not update fields"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    // Force linker to include zoltan symbols. This section is only needed since
    // Zoltan is a static library
    #ifdef FOAM_USE_ZOLTAN
        Info<< "renumberMesh built with zoltan support." << nl << endl;
        (void)zoltanRenumber::typeName;
    #endif

    // Get times list
    instantList Times = runTime.times();

    // Set startTime and endTime depending on -time and -latestTime options
    #include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    #include "createNamedMesh.H"
    const word oldInstance = mesh.pointsInstance();

    const bool readDict = args.optionFound("dict");
    const bool doFrontWidth = args.optionFound("frontWidth");
    const bool overwrite = args.optionFound("overwrite");
    const bool fields = !args.optionFound("noFields");

    label band;
    scalar profile;
    scalar sumSqrIntersect;
    getBand
    (
        doFrontWidth,
        mesh.nCells(),
        mesh.faceOwner(),
        mesh.faceNeighbour(),
        band,
        profile,
        sumSqrIntersect
    );

    reduce(band, maxOp<label>());
    reduce(profile, sumOp<scalar>());
    scalar rmsFrontwidth = Foam::sqrt
    (
        returnReduce
        (
            sumSqrIntersect,
            sumOp<scalar>()
        )/mesh.globalData().nTotalCells()
    );

    Info<< "Mesh size: " << mesh.globalData().nTotalCells() << nl
        << "Before renumbering :" << nl
        << "    band           : " << band << nl
        << "    profile        : " << profile << nl;

    if (doFrontWidth)
    {
        Info<< "    rms frontwidth : " << rmsFrontwidth << nl;
    }

    Info<< endl;

    bool sortCoupledFaceCells = false;
    bool writeMaps = false;
    bool orderPoints = false;
    label blockSize = 0;
    bool renumberSets = true;

    // Construct renumberMethod
    autoPtr<dictionary> renumberDictPtr;
    autoPtr<renumberMethod> renumberPtr;

    if (readDict)
    {
        renumberDictPtr.reset
        (
            new dictionary(systemDict("renumberMeshDict", args, mesh))
        );
        const dictionary& renumberDict = renumberDictPtr();

        renumberPtr = renumberMethod::New(renumberDict);

        sortCoupledFaceCells = renumberDict.lookupOrDefault
        (
            "sortCoupledFaceCells",
            false
        );
        if (sortCoupledFaceCells)
        {
            Info<< "Sorting cells on coupled boundaries to be last." << nl
                << endl;
        }

        blockSize = renumberDict.lookupOrDefault("blockSize", 0);
        if (blockSize > 0)
        {
            Info<< "Ordering cells into regions of size " << blockSize
                << " (using decomposition);"
                << " ordering faces into region-internal and region-external."
                << nl << endl;

            if (blockSize < 0 || blockSize >= mesh.nCells())
            {
                FatalErrorInFunction
                    << "Block size " << blockSize
                    << " should be positive integer"
                    << " and less than the number of cells in the mesh."
                    << exit(FatalError);
            }
        }

        orderPoints = renumberDict.lookupOrDefault("orderPoints", false);
        if (orderPoints)
        {
            Info<< "Ordering points into internal and boundary points." << nl
                << endl;
        }

        renumberDict.lookup("writeMaps") >> writeMaps;
        if (writeMaps)
        {
            Info<< "Writing renumber maps (new to old) to polyMesh." << nl
                << endl;
        }

        renumberSets = renumberDict.lookupOrDefault("renumberSets", true);
    }
    else
    {
        Info<< "Using default renumberMethod." << nl << endl;
        dictionary renumberDict;
        renumberPtr.reset(new CuthillMcKeeRenumber(renumberDict));
    }

    Info<< "Selecting renumberMethod " << renumberPtr().type() << nl << endl;



    // Read parallel reconstruct maps
    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );
    labelIOList pointProcAddressing
    (
        IOobject
        (
            "pointProcAddressing",
            mesh.pointsInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );
    labelIOList boundaryProcAddressing
    (
        IOobject
        (
            "boundaryProcAddressing",
            mesh.pointsInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        labelList(0)
    );


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    if (fields) Info<< "Reading geometric fields" << nl << endl;

    #include "readVolFields.H"
    #include "readSurfaceFields.H"
    #include "readPointFields.H"

    Info<< endl;

    // From renumbering:
    // - from new cell/face back to original cell/face
    labelList cellOrder;
    labelList faceOrder;

    if (blockSize > 0)
    {
        // Renumbering in two phases. Should be done in one so mapping of
        // fields is done correctly!

        label nBlocks = mesh.nCells()/blockSize;
        Info<< "nBlocks   = " << nBlocks << endl;

        // Read decompositionMethod dictionary
        dictionary decomposeDict(renumberDictPtr().subDict("blockCoeffs"));
        decomposeDict.set("numberOfSubdomains", nBlocks);

        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        autoPtr<decompositionMethod> decomposePtr = decompositionMethod::New
        (
            decomposeDict
        );

        labelList cellToRegion
        (
            decomposePtr().decompose
            (
                mesh,
                mesh.cellCentres()
            )
        );

        // Restore state
        UPstream::parRun() = oldParRun;

        // For debugging: write out region
        createScalarField
        (
            mesh,
            "cellDist",
            cellToRegion
        )().write();

        Info<< nl << "Written decomposition as volScalarField to "
            << "cellDist for use in postprocessing."
            << nl << endl;


        cellOrder = regionRenumber(renumberPtr(), mesh, cellToRegion);

        // Determine new to old face order with new cell numbering
        faceOrder = getRegionFaceOrder
        (
            mesh,
            cellOrder,
            cellToRegion
        );
    }
    else
    {
        // Determines sorted back to original cell ordering
        cellOrder = renumberPtr().renumber
        (
            mesh,
            mesh.cellCentres()
        );

        if (sortCoupledFaceCells)
        {
            // Change order so all coupled patch faceCells are at the end.
            const polyBoundaryMesh& pbm = mesh.boundaryMesh();

            // Collect all boundary cells on coupled patches
            label nBndCells = 0;
            forAll(pbm, patchi)
            {
                if (pbm[patchi].coupled())
                {
                    nBndCells += pbm[patchi].size();
                }
            }

            labelList reverseCellOrder = invert(mesh.nCells(), cellOrder);

            labelList bndCells(nBndCells);
            labelList bndCellMap(nBndCells);
            nBndCells = 0;
            forAll(pbm, patchi)
            {
                if (pbm[patchi].coupled())
                {
                    const labelUList& faceCells = pbm[patchi].faceCells();
                    forAll(faceCells, i)
                    {
                        label celli = faceCells[i];

                        if (reverseCellOrder[celli] != -1)
                        {
                            bndCells[nBndCells] = celli;
                            bndCellMap[nBndCells++] = reverseCellOrder[celli];
                            reverseCellOrder[celli] = -1;
                        }
                    }
                }
            }
            bndCells.setSize(nBndCells);
            bndCellMap.setSize(nBndCells);

            // Sort
            labelList order;
            sortedOrder(bndCellMap, order);

            // Redo newReverseCellOrder
            labelList newReverseCellOrder(mesh.nCells(), -1);

            label sortedI = mesh.nCells();
            forAllReverse(order, i)
            {
                label origCelli = bndCells[order[i]];
                newReverseCellOrder[origCelli] = --sortedI;
            }

            Info<< "Ordered all " << nBndCells << " cells with a coupled face"
                << " to the end of the cell list, starting at " << sortedI
                << endl;

            // Compact
            sortedI = 0;
            forAll(cellOrder, newCelli)
            {
                label origCelli = cellOrder[newCelli];
                if (newReverseCellOrder[origCelli] == -1)
                {
                    newReverseCellOrder[origCelli] = sortedI++;
                }
            }

            // Update sorted back to original (unsorted) map
            cellOrder = invert(mesh.nCells(), newReverseCellOrder);
        }


        // Determine new to old face order with new cell numbering
        faceOrder = getFaceOrder
        (
            mesh,
            cellOrder      // New to old cell
        );
    }


    if (!overwrite)
    {
        runTime++;
    }


    // Change the mesh.
    autoPtr<mapPolyMesh> map = reorderMesh(mesh, cellOrder, faceOrder);


    if (orderPoints)
    {
        polyTopoChange meshMod(mesh);
        autoPtr<mapPolyMesh> pointOrderMap = meshMod.changeMesh
        (
            mesh,
            false,      // inflate
            true,       // syncParallel
            false,      // orderCells
            orderPoints // orderPoints
        );

        // Combine point reordering into map.
        const_cast<labelList&>(map().pointMap()) = UIndirectList<label>
        (
            map().pointMap(),
            pointOrderMap().pointMap()
        )();

        inplaceRenumber
        (
            pointOrderMap().reversePointMap(),
            const_cast<labelList&>(map().reversePointMap())
        );
    }


    // Update fields
    mesh.updateMesh(map);

    // Update proc maps
    if
    (
        cellProcAddressing.headerOk()
     && cellProcAddressing.size() == mesh.nCells()
    )
    {
        Info<< "Renumbering processor cell decomposition map "
            << cellProcAddressing.name() << endl;

        cellProcAddressing = labelList
        (
            UIndirectList<label>(cellProcAddressing, map().cellMap())
        );
    }
    if
    (
        faceProcAddressing.headerOk()
     && faceProcAddressing.size() == mesh.nFaces()
    )
    {
        Info<< "Renumbering processor face decomposition map "
            << faceProcAddressing.name() << endl;

        faceProcAddressing = labelList
        (
            UIndirectList<label>(faceProcAddressing, map().faceMap())
        );

        // Detect any flips.
        const labelHashSet& fff = map().flipFaceFlux();
        forAllConstIter(labelHashSet, fff, iter)
        {
            label facei = iter.key();
            label masterFacei = faceProcAddressing[facei];

            faceProcAddressing[facei] = -masterFacei;

            if (masterFacei == 0)
            {
                FatalErrorInFunction
                    << " masterFacei:" << masterFacei << exit(FatalError);
            }
        }
    }
    if
    (
        pointProcAddressing.headerOk()
     && pointProcAddressing.size() == mesh.nPoints()
    )
    {
        Info<< "Renumbering processor point decomposition map "
            << pointProcAddressing.name() << endl;

        pointProcAddressing = labelList
        (
            UIndirectList<label>(pointProcAddressing, map().pointMap())
        );
    }


    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }


    {
        label band;
        scalar profile;
        scalar sumSqrIntersect;
        getBand
        (
            doFrontWidth,
            mesh.nCells(),
            mesh.faceOwner(),
            mesh.faceNeighbour(),
            band,
            profile,
            sumSqrIntersect
        );
        reduce(band, maxOp<label>());
        reduce(profile, sumOp<scalar>());
        scalar rmsFrontwidth = Foam::sqrt
        (
            returnReduce
            (
                sumSqrIntersect,
                sumOp<scalar>()
            )/mesh.globalData().nTotalCells()
        );

        Info<< "After renumbering :" << nl
            << "    band           : " << band << nl
            << "    profile        : " << profile << nl;

        if (doFrontWidth)
        {

            Info<< "    rms frontwidth : " << rmsFrontwidth << nl;
        }

        Info<< endl;
    }

    if (orderPoints)
    {
        // Force edge calculation (since only reason that points would need to
        // be sorted)
        (void)mesh.edges();

        label nTotPoints = returnReduce
        (
            mesh.nPoints(),
            sumOp<label>()
        );
        label nTotIntPoints = returnReduce
        (
            mesh.nInternalPoints(),
            sumOp<label>()
        );

        label nTotEdges = returnReduce
        (
            mesh.nEdges(),
            sumOp<label>()
        );
        label nTotIntEdges = returnReduce
        (
            mesh.nInternalEdges(),
            sumOp<label>()
        );
        label nTotInt0Edges = returnReduce
        (
            mesh.nInternal0Edges(),
            sumOp<label>()
        );
        label nTotInt1Edges = returnReduce
        (
            mesh.nInternal1Edges(),
            sumOp<label>()
        );

        Info<< "Points:" << nl
            << "    total   : " << nTotPoints << nl
            << "    internal: " << nTotIntPoints << nl
            << "    boundary: " << nTotPoints-nTotIntPoints << nl
            << "Edges:" << nl
            << "    total   : " << nTotEdges << nl
            << "    internal: " << nTotIntEdges << nl
            << "        internal using 0 boundary points: "
            << nTotInt0Edges << nl
            << "        internal using 1 boundary points: "
            << nTotInt1Edges-nTotInt0Edges << nl
            << "        internal using 2 boundary points: "
            << nTotIntEdges-nTotInt1Edges << nl
            << "    boundary: " << nTotEdges-nTotIntEdges << nl
            << endl;
    }


    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh to " << mesh.facesInstance() << endl;

    mesh.write();

    if (cellProcAddressing.headerOk())
    {
        cellProcAddressing.instance() = mesh.facesInstance();
        if (cellProcAddressing.size() == mesh.nCells())
        {
            cellProcAddressing.write();
        }
        else
        {
            // procAddressing file no longer valid. Delete it.
            const fileName fName(cellProcAddressing.filePath());
            if (fName.size())
            {
                Info<< "Deleting inconsistent processor cell decomposition"
                    << " map " << fName << endl;
                rm(fName);
            }
        }
    }

    if (faceProcAddressing.headerOk())
    {
        faceProcAddressing.instance() = mesh.facesInstance();
        if (faceProcAddressing.size() == mesh.nFaces())
        {
            faceProcAddressing.write();
        }
        else
        {
            const fileName fName(faceProcAddressing.filePath());
            if (fName.size())
            {
                Info<< "Deleting inconsistent processor face decomposition"
                    << " map " << fName << endl;
                rm(fName);
            }
        }
    }

    if (pointProcAddressing.headerOk())
    {
        pointProcAddressing.instance() = mesh.facesInstance();
        if (pointProcAddressing.size() == mesh.nPoints())
        {
            pointProcAddressing.write();
        }
        else
        {
            const fileName fName(pointProcAddressing.filePath());
            if (fName.size())
            {
                Info<< "Deleting inconsistent processor point decomposition"
                    << " map " << fName << endl;
                rm(fName);
            }
        }
    }

    if (boundaryProcAddressing.headerOk())
    {
        boundaryProcAddressing.instance() = mesh.facesInstance();
        if (boundaryProcAddressing.size() == mesh.boundaryMesh().size())
        {
            boundaryProcAddressing.write();
        }
        else
        {
            const fileName fName(boundaryProcAddressing.filePath());
            if (fName.size())
            {
                Info<< "Deleting inconsistent processor patch decomposition"
                    << " map " << fName << endl;
                rm(fName);
            }
        }
    }

    if (writeMaps)
    {
        // For debugging: write out region
        createScalarField
        (
            mesh,
            "origCellID",
            map().cellMap()
        )().write();

        createScalarField
        (
            mesh,
            "cellID",
            identity(mesh.nCells())
        )().write();

        Info<< nl << "Written current cellID and origCellID as volScalarField"
            << " for use in postprocessing."
            << nl << endl;

        labelIOList
        (
            IOobject
            (
                "cellMap",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            map().cellMap()
        ).write();

        labelIOList
        (
            IOobject
            (
                "faceMap",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            map().faceMap()
        ).write();

        labelIOList
        (
            IOobject
            (
                "pointMap",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            map().pointMap()
        ).write();
    }


    // Renumber sets if required
    if (renumberSets)
    {
        Info<< endl;

        // Read sets
        IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");

        {
            IOobjectList cSets(objects.lookupClass(cellSet::typeName));
            if (cSets.size())
            {
                Info<< "Renumbering cellSets:" << endl;
                forAllConstIter(IOobjectList, cSets, iter)
                {
                    cellSet cs(*iter());
                    Info<< "    " << cs.name() << endl;
                    cs.updateMesh(map());
                    cs.instance() = mesh.facesInstance();
                    cs.write();
                }
            }
        }

        {
            IOobjectList fSets(objects.lookupClass(faceSet::typeName));
            if (fSets.size())
            {
                Info<< "Renumbering faceSets:" << endl;
                forAllConstIter(IOobjectList, fSets, iter)
                {
                    faceSet fs(*iter());
                    Info<< "    " << fs.name() << endl;
                    fs.updateMesh(map());
                    fs.instance() = mesh.facesInstance();
                    fs.write();
                }
            }
        }

        {
            IOobjectList pSets(objects.lookupClass(pointSet::typeName));
            if (pSets.size())
            {
                Info<< "Renumbering pointSets:" << endl;
                forAllConstIter(IOobjectList, pSets, iter)
                {
                    pointSet ps(*iter());
                    Info<< "    " << ps.name() << endl;
                    ps.updateMesh(map());
                    ps.instance() = mesh.facesInstance();
                    ps.write();
                }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
