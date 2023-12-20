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

\*---------------------------------------------------------------------------*/

#include "domainDecomposition.H"
#include "fvMeshAdder.H"
#include "processorPolyPatch.H"
#include "processorCyclicPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::faceCoupleInfo>
Foam::domainDecomposition::determineCoupledFaces
(
    const label masterMeshProcStart,
    const label masterMeshProcEnd,
    const polyMesh& masterMesh,
    const label meshToAddProcStart,
    const label meshToAddProcEnd,
    const polyMesh& meshToAdd
)
{
    const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();
    const polyBoundaryMesh& addPatches = meshToAdd.boundaryMesh();

    DynamicList<label> masterFaces
    (
        masterMesh.nFaces() - masterMesh.nInternalFaces()
    );
    DynamicList<label> addFaces
    (
        meshToAdd.nFaces() - meshToAdd.nInternalFaces()
    );

    for
    (
        label masterProci = masterMeshProcStart;
        masterProci < masterMeshProcEnd;
        masterProci++
    )
    {
        for
        (
            label addProci = meshToAddProcStart;
            addProci < meshToAddProcEnd;
            addProci++
        )
        {
            const word masterToAddName
            (
                "procBoundary" + name(masterProci) + "to" + name(addProci)
            );
            const word addToMasterName
            (
                "procBoundary" + name(addProci) + "to" + name(masterProci)
            );

            const label masterToAddID =
                masterPatches.findIndex(masterToAddName);
            const label addToMasterID =
                addPatches.findIndex(addToMasterName);

            if (masterToAddID != -1 && addToMasterID != -1)
            {
                const polyPatch& masterPp = masterPatches[masterToAddID];

                forAll(masterPp, i)
                {
                    masterFaces.append(masterPp.start() + i);
                }

                const polyPatch& addPp = addPatches[addToMasterID];

                forAll(addPp, i)
                {
                    addFaces.append(addPp.start() + i);
                }
            }

            if ((masterToAddID != -1) != (addToMasterID != -1))
            {
                const label foundProci =
                    masterToAddID != -1 ? masterProci : addProci;
                const word& foundName =
                    masterToAddID != -1 ? masterToAddName : addToMasterName;

                const label missingProci =
                    masterToAddID != -1 ? addProci : masterProci;
                const word& missingName =
                    masterToAddID != -1 ? addToMasterName : masterToAddName;

                FatalErrorInFunction
                    << "Patch " << foundName << " found on processor "
                    << foundProci << " but corresponding patch "
                    << missingName << " missing on processor "
                    << missingProci << exit(FatalError);
            }
        }
    }

    masterFaces.shrink();
    addFaces.shrink();

    return autoPtr<faceCoupleInfo>
    (
        new faceCoupleInfo
        (
            masterMesh,
            masterFaces,
            meshToAdd,
            addFaces
        )
    );
}


void Foam::domainDecomposition::reconstruct()
{
    Info<< "Reconstructing meshes" << incrIndent << nl << endl;

    // ???
    PtrList<fvMesh> masterMeshes(nProcs());

    // ???
    for (label proci=0; proci<nProcs(); proci++)
    {
        masterMeshes.set
        (
            proci,
            new fvMesh
            (
                IOobject
                (
                    regionName_,
                    procMeshes()[0].facesInstance(),
                    runTimes_.completeTime()
                ),
                pointField(),
                faceList(),
                cellList()
            )
        );
        fvMesh& masterMesh = masterMeshes[proci];
        masterMesh.setPointsInstance(procMeshes()[0].pointsInstance());

        // Initialise the addressing
        procPointAddressing_[proci] =
            identityMap(procMeshes()[proci].nPoints());
        procFaceAddressing_[proci] = identityMap(procMeshes()[proci].nFaces());
        procCellAddressing_[proci] = identityMap(procMeshes()[proci].nCells());

        // Find shared points and faces
        autoPtr<faceCoupleInfo> couples = determineCoupledFaces
        (
            proci,
            proci,
            masterMesh,
            proci,
            proci,
            procMeshes()[proci]
        );

        // Add elements to the master mesh
        autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
        (
            masterMesh,
            procMeshes()[proci],
            couples,
            false
        );

        // Renumber addressing following the add
        inplaceRenumber
        (
            map().addedPointMap(),
            procPointAddressing_[proci]
        );
        inplaceRenumber
        (
            map().addedFaceMap(),
            procFaceAddressing_[proci]
        );
        inplaceRenumber
        (
            map().addedCellMap(),
            procCellAddressing_[proci]
        );
    }

    // Merge the meshes
    for (label step=2; step<nProcs()*2; step*=2)
    {
        for (label proci=0; proci<nProcs(); proci+=step)
        {
            label procj = proci + step/2;

            if (procj >= nProcs()) continue;

            Info<< indent << "Merging mesh " << proci
                << " with " << procj << endl;

            // Find shared points and faces
            autoPtr<faceCoupleInfo> couples = determineCoupledFaces
            (
                proci,
                procj,
                masterMeshes[proci],
                procj,
                proci+step,
                masterMeshes[procj]
            );

            // Add elements to mesh
            autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
            (
                masterMeshes[proci],
                masterMeshes[procj],
                couples,
                false
            );

            // Renumber processors that were already present
            for
            (
                label mergedProci = proci;
                mergedProci < procj;
                mergedProci ++
            )
            {
                inplaceRenumber
                (
                    map().oldPointMap(),
                    procPointAddressing_[mergedProci]
                );
                inplaceRenumber
                (
                    map().oldFaceMap(),
                    procFaceAddressing_[mergedProci]
                );
                inplaceRenumber
                (
                    map().oldCellMap(),
                    procCellAddressing_[mergedProci]
                );
            }

            // Renumber processors that were just added
            for
            (
                label addedProci = procj;
                addedProci < min(proci + step, nProcs());
                addedProci ++
            )
            {
                inplaceRenumber
                (
                    map().addedPointMap(),
                    procPointAddressing_[addedProci]
                );
                inplaceRenumber
                (
                    map().addedFaceMap(),
                    procFaceAddressing_[addedProci]
                );
                inplaceRenumber
                (
                    map().addedCellMap(),
                    procCellAddressing_[addedProci]
                );
            }

            masterMeshes.set(procj, nullptr);
        }
    }

    const polyBoundaryMesh& patches = masterMeshes[0].boundaryMesh();

    // Move all faces of processor cyclic patches into the associated cyclics
    if (!patches.findIndices<processorCyclicPolyPatch>().empty())
    {
        // Determine what patches are to be combined into each patch
        List<DynamicList<label>> patchPatches
        (
            patches.size(),
            DynamicList<label>(1)
        );
        forAll(patches, patchi)
        {
            patchPatches[patchi].append(patchi);
        }

        HashTable
        <
            label,
            FixedList<label, 3>,
            FixedList<label, 3>::Hash<>
        > nbrProcessorCyclicIDs;

        forAll(patches, patchi)
        {
            if (!isA<processorCyclicPolyPatch>(patches[patchi])) continue;

            const processorCyclicPolyPatch& pcpp =
                refCast<const processorCyclicPolyPatch>(patches[patchi]);

            const label proci = pcpp.myProcNo();
            const label nbrProci = pcpp.neighbProcNo();
            const label refPatchi = pcpp.referPatchIndex();
            const label nbrRefPatchi = pcpp.referPatch().nbrPatchIndex();

            FixedList<label, 3> key({proci, nbrProci, refPatchi});
            FixedList<label, 3> nbrKey({nbrProci, proci, nbrRefPatchi});

            if (nbrProcessorCyclicIDs.found(nbrKey))
            {
                const label nbrPatchi = nbrProcessorCyclicIDs[nbrKey];

                patchPatches[refPatchi].append
                (
                    patchPatches[patchi].remove()
                );
                patchPatches[nbrRefPatchi].append
                (
                    patchPatches[nbrPatchi].remove()
                );

                nbrProcessorCyclicIDs.erase(nbrKey);
            }
            else
            {
                nbrProcessorCyclicIDs.insert(key, patchi);
            }
        }

        // Build the new patch indices and a permutation map for the mesh faces
        labelList newPatchSizes(patches.size());
        labelList newPatchStarts(patches.size());
        labelList newToOldFace(masterMeshes[0].nFaces());
        SubList<label>(newToOldFace, masterMeshes[0].nInternalFaces()) =
            identityMap(masterMeshes[0].nInternalFaces());
        forAll(patches, patchi)
        {
            newPatchSizes[patchi] = 0;
            newPatchStarts[patchi] =
                patchi == 0
              ? masterMeshes[0].nInternalFaces()
              : newPatchStarts[patchi-1] + newPatchSizes[patchi-1];

            forAll(patchPatches[patchi], i)
            {
                const polyPatch& pp = patches[patchPatches[patchi][i]];

                SubList<label>
                (
                    newToOldFace,
                    pp.size(),
                    newPatchStarts[patchi] + newPatchSizes[patchi]
                ) = pp.start() + identityMap(pp.size());

                newPatchSizes[patchi] += pp.size();
            }
        }

        // Check that the face ordering is a permutation
        if (debug)
        {
            boolList newHasOldFace(newToOldFace.size(), false);
            forAll(newToOldFace, newFacei)
            {
                if (newHasOldFace[newToOldFace[newFacei]])
                {
                    FatalErrorInFunction
                        << "Face re-ordering is not valid"
                        << exit(FatalError);
                }
                newHasOldFace[newToOldFace[newFacei]] = true;
            }
        }

        // Modify the mesh
        const labelList oldToNewFace(invert(newToOldFace.size(), newToOldFace));
        masterMeshes[0].resetPrimitives
        (
            NullObjectMove<pointField>(),
            reorder(oldToNewFace, masterMeshes[0].faces()),
            reorder(oldToNewFace, masterMeshes[0].faceOwner()),
            reorder(oldToNewFace, masterMeshes[0].faceNeighbour()),
            newPatchSizes,
            newPatchStarts,
            true
        );

        // Update the addressing
        forAll(procFaceAddressing_, proci)
        {
            inplaceRenumber
            (
                oldToNewFace,
                procFaceAddressing_[proci]
            );
        }
    }

    // Filter out all processor patches
    {
        label nNonProcPatches = 0;
        labelList nonProcPatchIds(patches.size(), -1);

        forAll(masterMeshes[0].boundary(), patchi)
        {
            if (isA<processorPolyPatch>(patches[patchi]))
            {
                if (patches[patchi].size() != 0)
                {
                    FatalErrorInFunction
                        << "Non-empty processor patch \""
                        << patches[patchi].name()
                        << "\" found in reconstructed mesh"
                        << exit(FatalError);
                }
            }
            else
            {
                nonProcPatchIds[nNonProcPatches++] = patchi;
            }
        }

        nonProcPatchIds.resize(nNonProcPatches);

        masterMeshes[0].reorderPatches(nonProcPatchIds, false);
    }

    // Add turning index to the face addressing
    for (label proci=0; proci<nProcs(); proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        forAll(procFaceAddressing_[proci], procFacei)
        {
            const label completeFacei = procFaceAddressing_[proci][procFacei];

            if
            (
               !procMesh.isInternalFace(procFacei)
             && masterMeshes[0].isInternalFace(completeFacei)
            )
            {
                // Processor face is external, but the corresponding complete
                // face is internal. Turn if the owner has changed.

                const label procOwnCelli =
                    procMesh.faceOwner()[procFacei];
                const label completeOwnCelli =
                    masterMeshes[0].faceOwner()[completeFacei];

                if
                (
                    procCellAddressing_[proci][procOwnCelli]
                 == completeOwnCelli
                )
                {
                    procFaceAddressing_[proci][procFacei] ++;
                }
                else
                {
                    procFaceAddressing_[proci][procFacei] =
                        -1 - procFaceAddressing_[proci][procFacei];
                }
            }
            else
            {
                // Processor face is the same (internal/external) as the
                // corresponding complete face

                procFaceAddressing_[proci][procFacei] ++;
            }
        }
    }

    // Clear (and thus trigger re-generation) of finite volume face addressing
    procFaceAddressingBf_.clear();

    // Construct complete cell to proc map
    cellProc_.resize(masterMeshes[0].nCells());
    forAll(procCellAddressing_, proci)
    {
        forAll(procCellAddressing_[proci], procCelli)
        {
            cellProc_[procCellAddressing_[proci][procCelli]] = proci;
        }
    }

    // Set the complete mesh
    completeMesh_.reset(masterMeshes.set(0, nullptr).ptr());
    completeMesh_->setInstance(procMeshes()[0].facesInstance());
    completeMesh_->setPointsInstance(procMeshes()[0].pointsInstance());

    Info<< decrIndent;
}


void Foam::domainDecomposition::reconstructPoints()
{
    const label pointsCompare =
        compareInstances
        (
            completeMesh().pointsInstance(),
            procMeshes_[0].pointsInstance()
        );

    if (pointsCompare == 1)
    {
        pointField completePoints(completeMesh().nPoints());

        for (label proci = 0; proci < nProcs(); proci++)
        {
            const fvMesh& procMesh = procMeshes_[proci];

            completePoints.rmap(procMesh.points(), procPointAddressing_[proci]);
        }

        completeMesh_->setPoints(completePoints);
    }
}


// ************************************************************************* //
