/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "domainDecomposition.H"
#include "decompositionMethod.H"
#include "fvFieldDecomposer.H"
#include "IOobjectList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "hexRef8Data.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(domainDecomposition, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::domainDecomposition::mark
(
    const labelList& zoneElems,
    const label zoneI,
    labelList& elementToZone
)
{
    forAll(zoneElems, i)
    {
        label pointi = zoneElems[i];

        if (elementToZone[pointi] == -1)
        {
            // First occurrence
            elementToZone[pointi] = zoneI;
        }
        else if (elementToZone[pointi] >= 0)
        {
            // Multiple zones
            elementToZone[pointi] = -2;
        }
    }
}


void Foam::domainDecomposition::addInterProcFace
(
    const label facei,
    const label ownerProc,
    const label nbrProc,
    List<Map<label>>& nbrToInterPatch,
    List<DynamicList<DynamicList<label>>>& interPatchFaces
) const
{
    Map<label>::iterator patchiter = nbrToInterPatch[ownerProc].find(nbrProc);

    // Introduce turning index only for internal faces (are duplicated).
    label ownerIndex = facei+1;
    label nbrIndex = -(facei+1);

    if (patchiter != nbrToInterPatch[ownerProc].end())
    {
        // Existing interproc patch. Add to both sides.
        label toNbrProcPatchi = patchiter();
        interPatchFaces[ownerProc][toNbrProcPatchi].append(ownerIndex);

        if (mesh_.isInternalFace(facei))
        {
            label toOwnerProcPatchi = nbrToInterPatch[nbrProc][ownerProc];
            interPatchFaces[nbrProc][toOwnerProcPatchi].append(nbrIndex);
        }
    }
    else
    {
        // Create new interproc patches.
        label toNbrProcPatchi = nbrToInterPatch[ownerProc].size();
        nbrToInterPatch[ownerProc].insert(nbrProc, toNbrProcPatchi);
        DynamicList<label> oneFace;
        oneFace.append(ownerIndex);
        interPatchFaces[ownerProc].append(oneFace);

        if (mesh_.isInternalFace(facei))
        {
            label toOwnerProcPatchi = nbrToInterPatch[nbrProc].size();
            nbrToInterPatch[nbrProc].insert(ownerProc, toOwnerProcPatchi);
            oneFace.clear();
            oneFace.append(nbrIndex);
            interPatchFaces[nbrProc].append(oneFace);
        }
    }
}


Foam::labelList Foam::domainDecomposition::distributeCells()
{
    Info<< "\nCalculating distribution of cells" << endl;

    cpuTime decompositionTime;

    const dictionary decomposeParDict
    (
        decompositionMethod::decomposeParDict(mesh_.time())
    );

    scalarField cellWeights;
    if (decomposeParDict.found("weightField"))
    {
        const word weightName = decomposeParDict.lookup("weightField");

        volScalarField weights
        (
            IOobject
            (
                weightName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );
        cellWeights = weights.primitiveField();
    }

    const labelList result =
        decompositionMethod::NewDecomposer(decomposeParDict)->decompose
        (
            mesh_,
            cellWeights
        );

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;

    return result;
}


template<class BinaryOp>
inline void Foam::domainDecomposition::processInterCyclics
(
    const labelList& cellToProc,
    const polyBoundaryMesh& patches,
    List<DynamicList<DynamicList<label>>>& interPatchFaces,
    List<Map<label>>& procNbrToInterPatch,
    List<labelListList>& subPatchIDs,
    List<labelListList>& subPatchStarts,
    bool owner,
    BinaryOp bop
) const
{
    // Processor boundaries from split cyclics
    forAll(patches, patchi)
    {
        if (isA<cyclicPolyPatch>(patches[patchi]))
        {
            const cyclicPolyPatch& pp = refCast<const cyclicPolyPatch>
            (
                patches[patchi]
            );

            if (pp.owner() != owner)
            {
                continue;
            }

            // cyclic: check opposite side on this processor
            const labelUList& patchFaceCells = pp.faceCells();
            const labelUList& nbrPatchFaceCells =
                pp.nbrPatch().faceCells();

            // Store old sizes. Used to detect which inter-proc patches
            // have been added to.
            labelListList oldInterfaceSizes(nProcs_);
            forAll(oldInterfaceSizes, proci)
            {
                labelList& curOldSizes = oldInterfaceSizes[proci];

                curOldSizes.setSize(interPatchFaces[proci].size());
                forAll(curOldSizes, interI)
                {
                    curOldSizes[interI] =
                        interPatchFaces[proci][interI].size();
                }
            }

            // Add faces with different owner and neighbour processors
            forAll(patchFaceCells, facei)
            {
                const label ownerProc = cellToProc[patchFaceCells[facei]];
                const label nbrProc = cellToProc[nbrPatchFaceCells[facei]];
                if (bop(ownerProc, nbrProc))
                {
                    // inter - processor patch face found.
                    addInterProcFace
                    (
                        pp.start()+facei,
                        ownerProc,
                        nbrProc,
                        procNbrToInterPatch,
                        interPatchFaces
                    );
                }
            }

            // 1. Check if any faces added to existing interfaces
            forAll(oldInterfaceSizes, proci)
            {
                const labelList& curOldSizes = oldInterfaceSizes[proci];

                forAll(curOldSizes, interI)
                {
                    label oldSz = curOldSizes[interI];
                    if (interPatchFaces[proci][interI].size() > oldSz)
                    {
                        // Added faces to this interface. Add an entry
                        subPatchIDs[proci][interI].append(patchi);
                        subPatchStarts[proci][interI].append(oldSz);
                    }
                }
            }

            // 2. Any new interfaces
            forAll(subPatchIDs, proci)
            {
                label nIntfcs = interPatchFaces[proci].size();
                subPatchIDs[proci].setSize(nIntfcs, labelList(1, patchi));
                subPatchStarts[proci].setSize(nIntfcs, labelList(1, label(0)));
            }
        }
    }
}


void Foam::domainDecomposition::validate() const
{
    if (!procMeshes_.set(0))
    {
        FatalErrorInFunction
            << "Decomposition data requested but decomposition has not been "
            << "generated or read" << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::domainDecomposition::domainDecomposition
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    nProcs_
    (
        decompositionMethod::decomposeParDict(mesh_.time())
       .lookup<int>("numberOfSubdomains")
    ),
    facesInstancePointsPtr_(nullptr),
    distributed_(false),
    procPointAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procCellAddressing_(nProcs_),
    procRunTimes_(nProcs_),
    procMeshes_(nProcs_)
{
    readPoints();

    decompositionMethod::decomposeParDict(mesh_.time()).readIfPresent
    (
        "distributed",
        distributed_
    );

    // Create root databases for the processors
    for (label proci = 0; proci < nProcs_; proci++)
    {
        procRunTimes_.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                mesh_.time().rootPath(),
                mesh_.time().caseName()
               /fileName(word("processor") + Foam::name(proci)),
                word("system"),
                word("constant")
            )
        );

        procRunTimes_[proci].setTime(mesh_.time());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::domainDecomposition::~domainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::domainDecomposition::decompose()
{
    // Decide which cell goes to which processor
    const labelList cellToProc = distributeCells();

    // Distribute the cells according to the given processor label

    // calculate the addressing information for the original mesh
    Info<< "\nCalculating original mesh data" << endl;

    // set references to the original mesh
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const faceList& fcs = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    // loop through the list of processor labels for the cell and add the
    // cell shape to the list of cells for the appropriate processor

    Info<< "\nDistributing cells to processors" << endl;

    // Cells per processor
    procCellAddressing_ = invertOneToMany(nProcs_, cellToProc);

    Info<< "\nDistributing faces to processors" << endl;

    // Loop through all internal faces and decide which processor they belong to
    // First visit all internal faces. If cells at both sides belong to the
    // same processor, the face is an internal face. If they are different,
    // it belongs to both processors.

    // Internal faces
    forAll(procFaceAddressing_, proci)
    {
        procFaceAddressing_[proci].clear();
    }
    forAll(neighbour, facei)
    {
        if (cellToProc[owner[facei]] == cellToProc[neighbour[facei]])
        {
            // Face internal to processor. Notice no turning index.
            procFaceAddressing_[cellToProc[owner[facei]]].append(facei+1);
        }
    }

    // for all processors, set the size of start index and patch size
    // lists to the number of patches in the mesh
    labelListList procPatchSize(nProcs_);
    labelListList procPatchStartIndex(nProcs_);
    forAll(procPatchSize, proci)
    {
        procPatchSize[proci].setSize(patches.size());
        procPatchStartIndex[proci].setSize(patches.size());
    }

    forAll(patches, patchi)
    {
        // Reset size and start index for all processors
        forAll(procPatchSize, proci)
        {
            procPatchSize[proci][patchi] = 0;
            procPatchStartIndex[proci][patchi] =
                procFaceAddressing_[proci].size();
        }

        const label patchStart = patches[patchi].start();

        if (!isA<cyclicPolyPatch>(patches[patchi]))
        {
            // Normal patch. Add faces to processor where the cell
            // next to the face lives

            const labelUList& patchFaceCells =
                patches[patchi].faceCells();

            forAll(patchFaceCells, facei)
            {
                const label curProc = cellToProc[patchFaceCells[facei]];

                // add the face without turning index
                procFaceAddressing_[curProc].append(patchStart+facei+1);

                // increment the number of faces for this patch
                procPatchSize[curProc][patchi]++;
            }
        }
        else
        {
            const cyclicPolyPatch& pp = refCast<const cyclicPolyPatch>
            (
                patches[patchi]
            );
            // cyclic: check opposite side on this processor
            const labelUList& patchFaceCells = pp.faceCells();

            const labelUList& nbrPatchFaceCells =
                pp.nbrPatch().faceCells();

            forAll(patchFaceCells, facei)
            {
                const label curProc = cellToProc[patchFaceCells[facei]];
                const label nbrProc = cellToProc[nbrPatchFaceCells[facei]];
                if (curProc == nbrProc)
                {
                    // add the face without turning index
                    procFaceAddressing_[curProc].append(patchStart+facei+1);
                    // increment the number of faces for this patch
                    procPatchSize[curProc][patchi]++;
                }
            }
        }
    }

    // Done internal bits of the new mesh and the ordinary patches.

    // Per processor, from neighbour processor to the inter-processor patch
    // that communicates with that neighbour
    List<Map<label>> procNbrToInterPatch(nProcs_);

    // Per processor the faces per inter-processor patch
    List<DynamicList<DynamicList<label>>> interPatchFaces(nProcs_);

    // Processor boundaries from internal faces
    forAll(neighbour, facei)
    {
        label ownerProc = cellToProc[owner[facei]];
        label nbrProc = cellToProc[neighbour[facei]];

        if (ownerProc != nbrProc)
        {
            // inter - processor patch face found.
            addInterProcFace
            (
                facei,
                ownerProc,
                nbrProc,

                procNbrToInterPatch,
                interPatchFaces
            );
        }
    }

    // Add the proper processor faces to the sub information. For faces
    // originating from internal faces this is always -1.
    List<labelListList> subPatchIDs(nProcs_);
    List<labelListList> subPatchStarts(nProcs_);
    forAll(interPatchFaces, proci)
    {
        label nInterfaces = interPatchFaces[proci].size();

        subPatchIDs[proci].setSize(nInterfaces, labelList(1, label(-1)));
        subPatchStarts[proci].setSize(nInterfaces, labelList(1, label(0)));
    }

    // Special handling needed for the case that multiple processor cyclic
    // patches are created on each local processor domain, e.g. if a 3x3 case
    // is decomposed using the decomposition:
    //
    //              | 1 | 0 | 2 |
    //  cyclic left | 2 | 0 | 1 | cyclic right
    //              | 2 | 0 | 1 |
    //
    // - processors 1 and 2 will both have pieces of both cyclic left- and
    //   right sub-patches present
    // - the interface patch faces are stored in a single list, where each
    //   sub-patch is referenced into the list using a patch start index and
    //   size
    // - if the patches are in order (in the boundary file) of left, right
    //   - processor 1 will send: left, right
    //   - processor 1 will need to receive in reverse order: right, left
    //   - similarly for processor 2
    // - the sub-patches are therefore generated in 4 passes of the patch lists
    //   1. add faces from owner patch where local proc i < nbr proc i
    //   2. add faces from nbr patch where local proc i < nbr proc i
    //   3. add faces from owner patch where local proc i > nbr proc i
    //   4. add faces from nbr patch where local proc i > nbr proc i

    processInterCyclics
    (
        cellToProc,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        true,
        lessOp<label>()
    );

    processInterCyclics
    (
        cellToProc,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        false,
        lessOp<label>()
    );

    processInterCyclics
    (
        cellToProc,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        false,
        greaterOp<label>()
    );

    processInterCyclics
    (
        cellToProc,
        patches,
        interPatchFaces,
        procNbrToInterPatch,
        subPatchIDs,
        subPatchStarts,
        true,
        greaterOp<label>()
    );

    // Sort inter-proc patch by neighbour
    labelListList procNeighbourProcessors(nProcs_);
    labelListList procProcessorPatchSize(nProcs_);
    labelListList procProcessorPatchStartIndex(nProcs_);
    List<labelListList> procProcessorPatchSubPatchIDs(nProcs_);
    List<labelListList> procProcessorPatchSubPatchStarts(nProcs_);
    labelList order;
    forAll(procNbrToInterPatch, proci)
    {
        label nInterfaces = procNbrToInterPatch[proci].size();

        procNeighbourProcessors[proci].setSize(nInterfaces);
        procProcessorPatchSize[proci].setSize(nInterfaces);
        procProcessorPatchStartIndex[proci].setSize(nInterfaces);
        procProcessorPatchSubPatchIDs[proci].setSize(nInterfaces);
        procProcessorPatchSubPatchStarts[proci].setSize(nInterfaces);

        // Get sorted neighbour processors
        const Map<label>& curNbrToInterPatch = procNbrToInterPatch[proci];
        labelList nbrs = curNbrToInterPatch.toc();

        sortedOrder(nbrs, order);

        DynamicList<DynamicList<label>>& curInterPatchFaces =
            interPatchFaces[proci];

        forAll(nbrs, i)
        {
            const label nbrProc = nbrs[i];
            const label interPatch = curNbrToInterPatch[nbrProc];

            procNeighbourProcessors[proci][i] = nbrProc;
            procProcessorPatchSize[proci][i] =
                curInterPatchFaces[interPatch].size();
            procProcessorPatchStartIndex[proci][i] =
                procFaceAddressing_[proci].size();

            // Add size as last element to substarts and transfer
            subPatchStarts[proci][interPatch].append
            (
                curInterPatchFaces[interPatch].size()
            );
            procProcessorPatchSubPatchIDs[proci][i].transfer
            (
                subPatchIDs[proci][interPatch]
            );
            procProcessorPatchSubPatchStarts[proci][i].transfer
            (
                subPatchStarts[proci][interPatch]
            );

            // And add all the face labels for interPatch
            DynamicList<label>& interPatchFaces =
                curInterPatchFaces[interPatch];

            forAll(interPatchFaces, j)
            {
                procFaceAddressing_[proci].append(interPatchFaces[j]);
            }

            interPatchFaces.clearStorage();
        }

        curInterPatchFaces.clearStorage();
        procFaceAddressing_[proci].shrink();
    }

    if (debug)
    {
        forAll(procPatchStartIndex, proci)
        {
            Info<< "Processor:" << proci << endl;
            Info<< "    total faces:" << procFaceAddressing_[proci].size()
                << endl;

            const labelList& curProcPatchStartIndex =
                procPatchStartIndex[proci];

            forAll(curProcPatchStartIndex, patchi)
            {
                Info<< "    patch:" << patchi
                    << "\tstart:" << curProcPatchStartIndex[patchi]
                    << "\tsize:" << procPatchSize[proci][patchi]
                    << endl;
            }
        }
        Info<< endl;

        forAll(procNeighbourProcessors, proci)
        {
            Info<< "Processor " << proci << endl;
            forAll(procNeighbourProcessors[proci], i)
            {
                Info<< "    nbr:" << procNeighbourProcessors[proci][i] << endl;
                Info<< "    size:" << procProcessorPatchSize[proci][i] << endl;
                Info<< "    start:" << procProcessorPatchStartIndex[proci][i]
                    << endl;
            }
        }
        Info<< endl;

        forAll(procFaceAddressing_, proci)
        {
            Info<< "Processor:" << proci << endl;
            Info<< "    faces:" << procFaceAddressing_[proci] << endl;
        }
    }

    Info<< "\nDistributing points to processors" << endl;

    // For every processor, loop through the list of faces for the processor.
    // For every face, loop through the list of points and mark the point as
    // used for the processor. Collect the list of used points for the
    // processor.

    forAll(procPointAddressing_, proci)
    {
        boolList pointLabels(mesh_.nPoints(), false);

        // Get reference to list of used faces
        const labelList& procFaceLabels = procFaceAddressing_[proci];

        forAll(procFaceLabels, facei)
        {
            // Because of the turning index, some labels may be negative
            const labelList& facePoints = fcs[mag(procFaceLabels[facei]) - 1];

            forAll(facePoints, pointi)
            {
                // Mark the point as used
                pointLabels[facePoints[pointi]] = true;
            }
        }

        // Collect the used points
        labelList& procPointLabels = procPointAddressing_[proci];

        procPointLabels.setSize(pointLabels.size());

        label nUsedPoints = 0;

        forAll(pointLabels, pointi)
        {
            if (pointLabels[pointi])
            {
                procPointLabels[nUsedPoints] = pointi;

                nUsedPoints++;
            }
        }

        // Reset the size of used points
        procPointLabels.setSize(nUsedPoints);
    }

    Info<< "\nConstructing processor meshes" << endl;

    // Mark point/faces/cells that are in zones.
    // -1   : not in zone
    // -2   : in multiple zones
    // >= 0 : in single given zone
    // This will give direct lookup of elements that are in a single zone
    // and we'll only have to revert back to searching through all zones
    // for the duplicate elements

    // Point zones
    labelList pointToZone(mesh_.points().size(), -1);
    forAll(mesh_.pointZones(), zoneI)
    {
        mark(mesh_.pointZones()[zoneI], zoneI, pointToZone);
    }

    // Face zones
    labelList faceToZone(mesh_.faces().size(), -1);
    forAll(mesh_.faceZones(), zoneI)
    {
        mark(mesh_.faceZones()[zoneI], zoneI, faceToZone);
    }

    // Cell zones
    labelList cellToZone(mesh_.nCells(), -1);
    forAll(mesh_.cellZones(), zoneI)
    {
        mark(mesh_.cellZones()[zoneI], zoneI, cellToZone);
    }

    // Initialise information for reporting
    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;

    // Generate the meshes
    for (label proci = 0; proci < nProcs_; proci++)
    {
        // Create processor points
        const labelList& curPointLabels = procPointAddressing_[proci];
        const pointField& meshPoints = mesh_.points();
        labelList pointLookup(mesh_.nPoints(), -1);
        pointField procPoints(curPointLabels.size());
        forAll(curPointLabels, pointi)
        {
            procPoints[pointi] = meshPoints[curPointLabels[pointi]];
            pointLookup[curPointLabels[pointi]] = pointi;
        }

        // Create processor faces
        const labelList& curFaceLabels = procFaceAddressing_[proci];
        const faceList& meshFaces = mesh_.faces();
        labelList faceLookup(mesh_.nFaces(), -1);
        faceList procFaces(curFaceLabels.size());
        forAll(curFaceLabels, facei)
        {
            // Mark the original face as used
            // Remember to decrement the index by one (turning index)
            label curF = mag(curFaceLabels[facei]) - 1;

            faceLookup[curF] = facei;

            // get the original face
            labelList origFaceLabels;

            if (curFaceLabels[facei] >= 0)
            {
                // face not turned
                origFaceLabels = meshFaces[curF];
            }
            else
            {
                origFaceLabels = meshFaces[curF].reverseFace();
            }

            // translate face labels into local point list
            face& procFaceLabels = procFaces[facei];

            procFaceLabels.setSize(origFaceLabels.size());

            forAll(origFaceLabels, pointi)
            {
                procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
            }
        }

        // Create processor cells
        const labelList& curCellLabels = procCellAddressing_[proci];
        const cellList& meshCells = mesh_.cells();
        cellList procCells(curCellLabels.size());
        forAll(curCellLabels, celli)
        {
            const labelList& origCellLabels = meshCells[curCellLabels[celli]];
            cell& curCell = procCells[celli];
            curCell.setSize(origCellLabels.size());
            forAll(origCellLabels, cellFacei)
            {
                curCell[cellFacei] = faceLookup[origCellLabels[cellFacei]];
            }
        }

        // Create a processor mesh without a boundary. Two situations:
        //
        // - Static and topo-changing cases. Points and faces come from the
        //   same time/instance. The mesh will get constructed in the same
        //   instance.
        //
        // - Moving mesh cases. Points and faces come from different times. We
        //   read the points belonging to the faces instance. We then proceed
        //   as normal. Only at write time will we additionally write the
        //   current points.
        //
        if (mesh_.pointsInstance() != mesh_.facesInstance())
        {
            // Construct mesh from facesInstance.
            pointField facesInstancePoints
            (
                facesInstancePointsPtr_(),
                curPointLabels
            );

            procMeshes_.set
            (
                proci,
                new fvMesh
                (
                    IOobject
                    (
                        mesh_.polyMesh::name(), // region of undecomposed mesh
                        mesh_.facesInstance(),
                        procRunTimes_[proci]
                    ),
                    move(facesInstancePoints),
                    move(procFaces),
                    move(procCells)
                )
            );
        }
        else
        {
            procMeshes_.set
            (
                proci,
                new fvMesh
                (
                    IOobject
                    (
                        mesh_.polyMesh::name(), // region of undecomposed mesh
                        mesh_.facesInstance(),
                        procRunTimes_[proci]
                    ),
                    move(procPoints),
                    move(procFaces),
                    move(procCells)
                )
            );
        }
        fvMesh& procMesh = procMeshes_[proci];

        // Create processor boundary patch information
        const labelList& curPatchSizes = procPatchSize[proci];
        const labelList& curPatchStarts = procPatchStartIndex[proci];
        const labelList& curNeighbourProcessors =
            procNeighbourProcessors[proci];
        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize[proci];
        const labelList& curProcessorPatchStarts =
            procProcessorPatchStartIndex[proci];
        const labelListList& curSubPatchIDs =
            procProcessorPatchSubPatchIDs[proci];
        const labelListList& curSubStarts =
            procProcessorPatchSubPatchStarts[proci];
        const polyPatchList& meshPatches = mesh_.boundaryMesh();

        // Count the number of inter-proc patches
        label nInterProcPatches = 0;
        forAll(curSubPatchIDs, procPatchi)
        {
            nInterProcPatches += curSubPatchIDs[procPatchi].size();
        }
        List<polyPatch*> procPatches
        (
            curPatchSizes.size() + nInterProcPatches,
            nullptr
        );

        label nPatches = 0;

        // Map existing non-proc patches
        forAll(curPatchSizes, patchi)
        {
            // Get the face labels consistent with the field mapping
            // (reuse the patch field mappers)
            const polyPatch& meshPatch = meshPatches[patchi];

            fvFieldDecomposer::patchFieldDecomposer patchMapper
            (
                SubList<label>
                (
                    curFaceLabels,
                    curPatchSizes[patchi],
                    curPatchStarts[patchi]
                ),
                meshPatch.start()
            );

            // Map existing patches
            procPatches[nPatches] = meshPatch.clone
            (
                procMesh.boundaryMesh(),
                nPatches,
                patchMapper.addressing(),
                curPatchStarts[patchi]
            ).ptr();

            nPatches++;
        }

        // Create new inter-proc patches
        forAll(curProcessorPatchSizes, procPatchi)
        {
            const labelList& subPatchID = curSubPatchIDs[procPatchi];
            const labelList& subStarts = curSubStarts[procPatchi];

            label curStart = curProcessorPatchStarts[procPatchi];

            forAll(subPatchID, i)
            {
                const label size =
                    i < subPatchID.size()-1
                  ? subStarts[i+1] - subStarts[i]
                  : curProcessorPatchSizes[procPatchi] - subStarts[i];

                if (subPatchID[i] == -1)
                {
                    // From internal faces
                    procPatches[nPatches] =
                        new processorPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi]
                        );
                }
                else
                {
                    const coupledPolyPatch& cpp =
                        refCast<const coupledPolyPatch>
                        (mesh_.boundaryMesh()[subPatchID[i]]);

                    procPatches[nPatches] =
                        new processorCyclicPolyPatch
                        (
                            size,
                            curStart,
                            nPatches,
                            procMesh.boundaryMesh(),
                            proci,
                            curNeighbourProcessors[procPatchi],
                            cpp.name()
                        );
                }

                curStart += size;

                nPatches++;
            }
        }

        // Add patches to the mesh
        procMesh.addFvPatches(procPatches);

        // Create point zones
        {
            const meshPointZones& pz = mesh_.pointZones();

            // Go through all the zoned points and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zonePoints(pz.size());

            // Estimate size
            forAll(zonePoints, zoneI)
            {
                zonePoints[zoneI].setCapacity(pz[zoneI].size() / nProcs_);
            }

            // Use the pointToZone map to find out the single zone (if any),
            // use slow search only for shared points.
            forAll(curPointLabels, pointi)
            {
                label curPoint = curPointLabels[pointi];

                label zoneI = pointToZone[curPoint];

                if (zoneI >= 0)
                {
                    // Single zone.
                    zonePoints[zoneI].append(pointi);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(pz, zoneI)
                    {
                        label index = pz[zoneI].whichPoint(curPoint);

                        if (index != -1)
                        {
                            zonePoints[zoneI].append(pointi);
                        }
                    }
                }
            }

            procMesh.pointZones().clearAddressing();
            procMesh.pointZones().setSize(zonePoints.size());
            forAll(zonePoints, zoneI)
            {
                procMesh.pointZones().set
                (
                    zoneI,
                    pz[zoneI].clone
                    (
                        procMesh.pointZones(),
                        zoneI,
                        zonePoints[zoneI].shrink()
                    )
                );
            }

            if (pz.size())
            {
                // Force writing on all processors
                procMesh.pointZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Create face zones
        {
            const meshFaceZones& fz = mesh_.faceZones();

            // Go through all the zoned face and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneFaces(fz.size());
            List<DynamicList<bool>> zoneFaceFlips(fz.size());

            // Estimate size
            forAll(zoneFaces, zoneI)
            {
                label procSize = fz[zoneI].size() / nProcs_;

                zoneFaces[zoneI].setCapacity(procSize);
                zoneFaceFlips[zoneI].setCapacity(procSize);
            }

            // Go through all the zoned faces and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll(curFaceLabels, facei)
            {
                // Remember to decrement the index by one (turning index)
                //
                label curF = mag(curFaceLabels[facei]) - 1;

                label zoneI = faceToZone[curF];

                if (zoneI >= 0)
                {
                    // Single zone. Add the face
                    zoneFaces[zoneI].append(facei);

                    label index = fz[zoneI].whichFace(curF);

                    bool flip = fz[zoneI].flipMap()[index];

                    if (curFaceLabels[facei] < 0)
                    {
                        flip = !flip;
                    }

                    zoneFaceFlips[zoneI].append(flip);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(fz, zoneI)
                    {
                        label index = fz[zoneI].whichFace(curF);

                        if (index != -1)
                        {
                            zoneFaces[zoneI].append(facei);

                            bool flip = fz[zoneI].flipMap()[index];

                            if (curFaceLabels[facei] < 0)
                            {
                                flip = !flip;
                            }

                            zoneFaceFlips[zoneI].append(flip);
                        }
                    }
                }
            }

            procMesh.faceZones().clearAddressing();
            procMesh.faceZones().setSize(zoneFaces.size());
            forAll(zoneFaces, zoneI)
            {
                procMesh.faceZones().set
                (
                    zoneI,
                    fz[zoneI].clone
                    (
                        zoneFaces[zoneI].shrink(),          // addressing
                        zoneFaceFlips[zoneI].shrink(),      // flipmap
                        zoneI,
                        procMesh.faceZones()
                    )
                );
            }

            if (fz.size())
            {
                // Force writing on all processors
                procMesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Create cell zones
        {
            const meshCellZones& cz = mesh_.cellZones();

            // Go through all the zoned cells and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            List<DynamicList<label>> zoneCells(cz.size());

            // Estimate size
            forAll(zoneCells, zoneI)
            {
                zoneCells[zoneI].setCapacity(cz[zoneI].size() / nProcs_);
            }

            forAll(curCellLabels, celli)
            {
                label curCelli = curCellLabels[celli];

                label zoneI = cellToZone[curCelli];

                if (zoneI >= 0)
                {
                    // Single zone.
                    zoneCells[zoneI].append(celli);
                }
                else if (zoneI == -2)
                {
                    // Multiple zones. Lookup.
                    forAll(cz, zoneI)
                    {
                        label index = cz[zoneI].whichCell(curCelli);

                        if (index != -1)
                        {
                            zoneCells[zoneI].append(celli);
                        }
                    }
                }
            }

            procMesh.cellZones().clearAddressing();
            procMesh.cellZones().setSize(zoneCells.size());
            forAll(zoneCells, zoneI)
            {
                procMesh.cellZones().set
                (
                    zoneI,
                    cz[zoneI].clone
                    (
                        zoneCells[zoneI].shrink(),
                        zoneI,
                        procMesh.cellZones()
                    )
                );
            }

            if (cz.size())
            {
                // Force writing on all processors
                procMesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
            }
        }

        // Report processor and update global statistics
        {
            Info<< endl
                << "Processor " << proci << nl
                << "    Number of cells = " << procMesh.nCells()
                << endl;

            maxProcCells = max(maxProcCells, procMesh.nCells());

            label nBoundaryFaces = 0;
            label nProcPatches = 0;
            label nProcFaces = 0;

            forAll(procMesh.boundaryMesh(), patchi)
            {
                if (isA<processorPolyPatch>(procMesh.boundaryMesh()[patchi]))
                {
                    const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>
                    (
                        procMesh.boundaryMesh()[patchi]
                    );

                    Info<< "    Number of faces shared with processor "
                        << ppp.neighbProcNo() << " = " << ppp.size() << endl;

                    nProcPatches++;
                    nProcFaces += ppp.size();
                }
                else
                {
                    nBoundaryFaces += procMesh.boundaryMesh()[patchi].size();
                }
            }

            Info<< "    Number of processor patches = " << nProcPatches << nl
                << "    Number of processor faces = " << nProcFaces << nl
                << "    Number of boundary faces = " << nBoundaryFaces << endl;

            totProcFaces += nProcFaces;
            totProcPatches += nProcPatches;
            maxProcPatches = max(maxProcPatches, nProcPatches);
            maxProcFaces = max(maxProcFaces, nProcFaces);
        }
    }

    // Determine the average number of processor elements
    scalar avgProcCells = scalar(mesh_.nCells())/nProcs_;
    scalar avgProcPatches = scalar(totProcPatches)/nProcs_;
    scalar avgProcFaces = scalar(totProcFaces)/nProcs_;

    // Prevent division by zero in the case of all faces on one processor
    if (totProcPatches == 0)
    {
        avgProcPatches = 1;
    }
    if (totProcFaces == 0)
    {
        avgProcFaces = 1;
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of cells = " << maxProcCells
        << " (" << 100.0*(maxProcCells-avgProcCells)/avgProcCells
        << "% above average " << avgProcCells << ")" << nl
        << "Max number of processor patches = " << maxProcPatches
        << " (" << 100.0*(maxProcPatches-avgProcPatches)/avgProcPatches
        << "% above average " << avgProcPatches << ")" << nl
        << "Max number of faces between processors = " << maxProcFaces
        << " (" << 100.0*(maxProcFaces-avgProcFaces)/avgProcFaces
        << "% above average " << avgProcFaces << ")" << nl
        << endl;
}


void Foam::domainDecomposition::setTime
(
    const instant& inst,
    const label newIndex
)
{
    for (label proci = 0; proci < nProcs_; proci++)
    {
        procRunTimes_[proci].setTime(inst, newIndex);
    }
}


void Foam::domainDecomposition::readPoints()
{
    if (mesh_.pointsInstance() != mesh_.facesInstance())
    {
        facesInstancePointsPtr_.reset
        (
            new pointIOField
            (
                IOobject
                (
                    "points",
                    mesh_.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
    }
}


void Foam::domainDecomposition::readAddressing()
{
    for (label proci = 0; proci < nProcs_; proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        procPointAddressing_[proci] =
            labelIOList
            (
                IOobject
                (
                    "pointProcAddressing",
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        procFaceAddressing_[proci] =
            labelIOList
            (
                IOobject
                (
                    "faceProcAddressing",
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        procCellAddressing_[proci] =
            labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
    }
}


void Foam::domainDecomposition::read()
{
    for (label proci = 0; proci < nProcs_; proci++)
    {
        procMeshes_.set
        (
            proci,
            new fvMesh
            (
                IOobject
                (
                    mesh_.polyMesh::name(), // region of undecomposed mesh
                    procRunTimes_[proci].timeName(),
                    procRunTimes_[proci]
                ),
                false
            )
        );
    }

    readAddressing();
}


Foam::fvMesh::readUpdateState Foam::domainDecomposition::readUpdate()
{
    fvMesh::readUpdateState stat = fvMesh::UNCHANGED;

    forAll(procRunTimes_, proci)
    {
        fvMesh::readUpdateState procStat = procMeshes_[proci].readUpdate();

        if (procStat > stat)
        {
            stat = procStat;
        }
    }

    if (mesh_.pointsInstance() != mesh_.facesInstance())
    {
        readPoints();
    }
    else
    {
        facesInstancePointsPtr_.clear();
    }

    if
    (
        stat == fvMesh::TOPO_CHANGE
     || stat == fvMesh::TOPO_PATCH_CHANGE
    )
    {
        readAddressing();
    }

    return stat;
}


Foam::labelList Foam::domainDecomposition::cellToProc() const
{
    labelList result(mesh_.nCells());

    forAll(procCellAddressing_, proci)
    {
        forAll(procCellAddressing_[proci], procCelli)
        {
            result[procCellAddressing_[proci][procCelli]] = proci;
        }
    }

    return result;
}


void Foam::domainDecomposition::writePoints() const
{
    for (label proci = 0; proci < nProcs_; proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        if (mesh_.pointsInstance() != mesh_.facesInstance())
        {
            pointField procPoints
            (
                mesh_.points(),
                procPointAddressing_[proci]
            );

            pointIOField pointsInstancePoints
            (
                IOobject
                (
                    "points",
                    mesh_.pointsInstance(),
                    polyMesh::meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                move(procPoints)
            );

            pointsInstancePoints.write();
        }
    }
}


void Foam::domainDecomposition::writeAddressing() const
{
    for (label proci = 0; proci < nProcs_; proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procPointAddressing_[proci]
        );
        pointProcAddressing.write();

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procFaceAddressing_[proci]
        );
        faceProcAddressing.write();

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procCellAddressing_[proci]
        );
        cellProcAddressing.write();
    }
}


void Foam::domainDecomposition::write(const bool decomposeSets) const
{
    validate();

    // Read sets
    PtrList<const cellSet> cellSets;
    PtrList<const faceSet> faceSets;
    PtrList<const pointSet> pointSets;
    if (decomposeSets)
    {
        IOobjectList objects(mesh_, mesh_.facesInstance(), "polyMesh/sets");
        {
            IOobjectList cSets(objects.lookupClass(cellSet::typeName));
            forAllConstIter(IOobjectList, cSets, iter)
            {
                cellSets.append(new cellSet(*iter()));
            }
        }
        {
            IOobjectList fSets(objects.lookupClass(faceSet::typeName));
            forAllConstIter(IOobjectList, fSets, iter)
            {
                faceSets.append(new faceSet(*iter()));
            }
        }
        {
            IOobjectList pSets(objects.lookupClass(pointSet::typeName));
            forAllConstIter(IOobjectList, pSets, iter)
            {
                pointSets.append(new pointSet(*iter()));
            }
        }
    }

    // Read refinement data (if any)
    hexRef8Data refinementData
    (
        IOobject
        (
            "dummy",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    // Write out the meshes
    for (label proci = 0; proci < nProcs_; proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        // Set the precision of the points data to be min 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        // Write the processor mesh
        procMesh.write();

        // Write any sets
        if (decomposeSets)
        {
            forAll(cellSets, i)
            {
                const cellSet& cs = cellSets[i];
                cellSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(procCellAddressing_[proci], i)
                {
                    if (cs.found(procCellAddressing_[proci][i]))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
            forAll(faceSets, i)
            {
                const faceSet& cs = faceSets[i];
                faceSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(procFaceAddressing_[proci], i)
                {
                    if (cs.found(mag(procFaceAddressing_[proci][i])-1))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
            forAll(pointSets, i)
            {
                const pointSet& cs = pointSets[i];
                pointSet set(procMesh, cs.name(), cs.size()/nProcs_);
                forAll(procPointAddressing_[proci], i)
                {
                    if (cs.found(procPointAddressing_[proci][i]))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
        }

        // Write refinement data (if any)
        hexRef8Data
        (
            IOobject
            (
                "dummy",
                mesh_.facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            refinementData,
            procCellAddressing_[proci],
            procPointAddressing_[proci]
        ).write();
    }

    // Write points if pointsInstance differing from facesInstance
    writePoints();

    // Write decomposition addressing
    writeAddressing();
}


// ************************************************************************* //
