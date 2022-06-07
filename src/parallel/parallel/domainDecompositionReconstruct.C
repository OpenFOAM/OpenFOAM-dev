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
#include "fvMeshAdder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

autoPtr<faceCoupleInfo> determineCoupledFaces
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
                masterPatches.findPatchID(masterToAddName);
            const label addToMasterID =
                addPatches.findPatchID(addToMasterName);

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

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
        procPointAddressing_[proci] = identity(procMeshes()[proci].nPoints());
        procFaceAddressing_[proci] = identity(procMeshes()[proci].nFaces());
        procCellAddressing_[proci] = identity(procMeshes()[proci].nCells());

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
            couples
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

            // Find shared points/faces
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
                couples
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

    // Set the complete mesh
    completeMesh_.reset(masterMeshes.set(0, nullptr).ptr());
    completeMesh_->setInstance(procMeshes()[0].facesInstance());
    completeMesh_->setPointsInstance(procMeshes()[0].pointsInstance());

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
             && completeMesh().isInternalFace(completeFacei)
            )
            {
                // Processor face is external, but the corresponding complete
                // face is internal. Turn if the owner has changed.

                const label procOwnCelli =
                    procMesh.faceOwner()[procFacei];
                const label completeOwnCelli =
                    completeMesh().faceOwner()[completeFacei];

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
    cellProc_.resize(completeMesh().nCells());
    forAll(procCellAddressing_, proci)
    {
        forAll(procCellAddressing_[proci], procCelli)
        {
            cellProc_[procCellAddressing_[proci][procCelli]] = proci;
        }
    }

    Info<< decrIndent << endl;
}


// ************************************************************************* //
