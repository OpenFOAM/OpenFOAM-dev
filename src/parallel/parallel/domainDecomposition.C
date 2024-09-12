/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "IOobjectList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "fvMeshStitcher.H"
#include "pointSet.H"
#include "hexRef8Data.H"
#include "cyclicFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "nonConformalFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(domainDecomposition, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::domainDecomposition::compareInstances
(
    const fileName& a,
    const fileName& b
) const
{
    const word& constant = runTimes_.completeTime().constant();

    if (a == constant && b == constant) return 0;

    if (a == constant) return +1;

    if (b == constant) return -1;

    const scalar aValue = instant(a).value();
    const scalar bValue = instant(b).value();

    if (aValue < bValue) return +1;

    if (aValue > bValue) return -1;

    return 0;
}


void Foam::domainDecomposition::validateComplete() const
{
    if (!completeMesh_.valid())
    {
        FatalErrorInFunction
            << "Complete data requested but complete mesh has not been "
            << "generated or read" << exit(FatalError);
    }
}


void Foam::domainDecomposition::validateProcs() const
{
    if (!procMeshes_.set(0))
    {
        FatalErrorInFunction
            << "Decomposed data requested but decomposed mesh has not been "
            << "generated or read" << exit(FatalError);
    }
}


void Foam::domainDecomposition::readComplete()
{
    completeMesh_.reset
    (
        new fvMesh
        (
            IOobject
            (
                regionName_,
                runTimes_.completeTime().name(),
                meshPath_,
                runTimes_.completeTime(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            false
        )
    );

    completeMesh_->postConstruct(false, fvMesh::stitchType::none);
}


void Foam::domainDecomposition::readProcs()
{
    for (label proci = 0; proci < nProcs(); proci++)
    {
        procMeshes_.set
        (
            proci,
            new fvMesh
            (
                IOobject
                (
                    regionName_,
                    runTimes_.procTimes()[proci].name(),
                    meshPath_,
                    runTimes_.procTimes()[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                false
            )
        );

        procMeshes_[proci].postConstruct(false, fvMesh::stitchType::none);
    }
}


void Foam::domainDecomposition::readCompleteAddressing()
{
    cellProc_ =
        labelIOList
        (
            IOobject
            (
                "cellProc",
                completeMesh().facesInstance(),
                completeMesh().meshSubDir,
                completeMesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
}


void Foam::domainDecomposition::readProcsAddressing()
{
    for (label proci = 0; proci < nProcs(); proci++)
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


void Foam::domainDecomposition::readAddressing()
{
    readCompleteAddressing();
    readProcsAddressing();
    procFaceAddressingBf_.clear();
}


Foam::fvMesh::readUpdateState Foam::domainDecomposition::readUpdate()
{
    validateComplete();
    validateProcs();

    // Do read-update on all meshes
    fvMesh::readUpdateState stat =
        completeMesh_->readUpdate(fvMesh::stitchType::none);

    forAll(runTimes_.procTimes(), proci)
    {
        fvMesh::readUpdateState procStat =
            procMeshes_[proci].readUpdate(fvMesh::stitchType::none);

        stat = procStat > stat ? procStat : stat;
    }

    return stat;
}


void Foam::domainDecomposition::writeCompleteAddressing() const
{
    labelIOList cellProc
    (
        IOobject
        (
            "cellProc",
            completeMesh().facesInstance(),
            completeMesh().meshSubDir,
            completeMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cellProc_
    );

    cellProc.write();
}


void Foam::domainDecomposition::writeProcsAddressing() const
{
    for (label proci = 0; proci < nProcs(); proci++)
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


void Foam::domainDecomposition::writeAddressing() const
{
    writeCompleteAddressing();
    writeProcsAddressing();
}


void Foam::domainDecomposition::writeProcPoints(const fileName& inst)
{
    IOobject completePointsIo
    (
        "points",
        inst,
        polyMesh::meshSubDir,
        completeMesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!completePointsIo.headerOk()) return;

    const pointIOField completePoints(completePointsIo);

    for (label proci = 0; proci < nProcs(); proci++)
    {
        pointIOField procPoints
        (
            IOobject
            (
                "points",
                inst,
                polyMesh::meshSubDir,
                procMeshes()[proci],
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointField
            (
                completePoints,
                procPointAddressing_[proci]
            )
        );

        procPoints.write();
    }
}


void Foam::domainDecomposition::writeCompletePoints(const fileName& inst)
{
    pointIOField completePoints
    (
        IOobject
        (
            "points",
            inst,
            polyMesh::meshSubDir,
            completeMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointField(completeMesh().nPoints())
    );

    for (label proci = 0; proci < nProcs(); proci++)
    {
        IOobject procPointsIo
        (
            "points",
            inst,
            polyMesh::meshSubDir,
            procMeshes()[proci],
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!procPointsIo.headerOk()) return;

        completePoints.rmap
        (
            pointIOField(procPointsIo),
            procPointAddressing_[proci]
        );
    }

    completePoints.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::domainDecomposition::domainDecomposition
(
    const processorRunTimes& runTimes,
    const fileName& meshPath,
    const word& regionName,
    const multiDomainDecomposition& regionMeshes
)
:
    runTimes_(runTimes),
    meshPath_(meshPath),
    regionName_(regionName),
    completeMesh_(nullptr),
    procMeshes_(nProcs()),
    regionMeshes_(regionMeshes),
    cellProc_(),
    procPointAddressing_(nProcs()),
    procFaceAddressing_(nProcs()),
    procCellAddressing_(nProcs()),
    procFaceAddressingBf_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::domainDecomposition::~domainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::domainDecomposition::readDecompose()
{
    readComplete();

    typeIOobject<labelIOList> addrIo
    (
        "cellProc",
        completeMesh().facesInstance(),
        polyMesh::meshSubDir,
        completeMesh()
    );
    IOobject procFaceIo
    (
        "faces",
        completeMesh().facesInstance(),
        completeMesh().meshDir(),
        runTimes_.procTimes()[0]
    );
    typeIOobject<labelIOList> procAddrIo
    (
        "cellProcAddressing",
        completeMesh().facesInstance(),
        completeMesh().meshDir(),
        runTimes_.procTimes()[0]
    );

    const bool load = addrIo.headerOk() && procFaceIo.headerOk();

    if (load)
    {
        readProcs();

        if (procAddrIo.headerOk())
        {
            readAddressing();
        }
        else
        {
            readCompleteAddressing();

            FatalErrorInFunction
                << nl << "    Processor meshes exist but have no addressing."
                << nl << nl << "    This could be because the processor meshes "
                << "have changed. Decomposing the" << nl << "    mesh would "
                << "overwrite that change. If you are sure that this is "
                << "appropriate," << nl << "    then delete the "
                << fileName("processor*")/procFaceIo.relativePath().c_str()
                << " directories and re-run this" << nl << "    command."
                << exit(FatalError);
        }

        decomposePoints();
    }
    else
    {
        if
        (
            completeMesh().facesInstance()
         != runTimes_.completeTime().name()
         && completeMesh().facesInstance()
         != runTimes_.completeTime().constant()
        )
        {
            FatalErrorInFunction
                << "Cannot begin mesh decomposition at time "
                << fileName(runTimes_.completeTime().name()) << nl
                << "The mesh at this instant is that of an earlier"
                << " time " << completeMesh().facesInstance() << nl
                << "Decomposition must start from this earlier time"
                << exit(FatalError);
        }

        decompose();
    }

    return !load;
}


void Foam::domainDecomposition::postReadDecompose()
{
    completeMesh_->stitcher().connect(false, false, true);
}


void Foam::domainDecomposition::unconformReadDecompose()
{
    if (!completeConformal())
    {
        procFaceAddressingBf_.clear();

        forAll(procMeshes_, proci)
        {
            procMeshes_[proci].conform();
        }

        unconform();
    }
}


void Foam::domainDecomposition::writeReadDecompose
(
    const bool decomposed,
    const bool doSets
)
{
    writeProcs(doSets);

    if (decomposed)
    {
        writeProcPoints(completeMesh().facesInstance());
    }
}


bool Foam::domainDecomposition::readReconstruct()
{
    readProcs();

    IOobject faceIo
    (
        "faces",
        procMeshes()[0].facesInstance(),
        procMeshes()[0].meshDir(),
        runTimes_.completeTime()
    );
    typeIOobject<labelIOList> addrIo
    (
        "cellProc",
        procMeshes()[0].facesInstance(),
        procMeshes()[0].meshDir(),
        runTimes_.completeTime()
    );
    typeIOobject<labelIOList> procAddrIo
    (
        "cellProcAddressing",
        procMeshes()[0].facesInstance(),
        polyMesh::meshSubDir,
        procMeshes()[0]
    );

    const bool load = faceIo.headerOk() && procAddrIo.headerOk();

    if (load)
    {
        typeIOobject<pointIOField> completePointsIo
        (
            "points",
            procMeshes()[0].pointsInstance(),
            procMeshes()[0].meshDir(),
            runTimes_.completeTime()
        );

        readComplete();

        if (addrIo.headerOk())
        {
            readAddressing();
        }
        else
        {
            readProcsAddressing();

            WarningInFunction
                << nl << "    A complete mesh exists but has no "
                << addrIo.name() << " addressing." << nl << nl << "    This "
                << "could be because the complete mesh has changed. "
                << "Reconstructing the" << nl << "    mesh would overwrite "
                << "that change. If you are sure that this is appropriate,"
                << nl << "    then delete the " << faceIo.relativePath()
                << " directory and re-run this command." << nl << nl
                << "    Or, it could be because the complete and processor "
                << "meshes were decomposed" << nl << "    by a version of "
                << "OpenFOAM that pre-dates the automatic generation of "
                << nl << "    " << addrIo.name() << " addressing. This will be "
                << "assumed and the " << addrIo.name() << " addressing will"
                << nl << "    be re-built" << nl << endl;

            cellProc_ = labelList(completeMesh().nCells(), -1);

            for (label proci = 0; proci < nProcs(); proci++)
            {
                UIndirectList<label>
                (
                    cellProc_,
                    procCellAddressing_[proci]
                ) = proci;
            }

            writeCompleteAddressing();
        }

        reconstructPoints();
    }
    else
    {
        if
        (
            procMeshes()[0].facesInstance()
         != runTimes_.procTimes()[0].name()
         && procMeshes()[0].facesInstance()
         != runTimes_.procTimes()[0].constant()
        )
        {
            FatalErrorInFunction
                << "Cannot begin mesh reconstruction at time "
                << fileName(runTimes_.procTimes()[0].name()) << nl
                << "The mesh at this instant is that of an earlier"
                << " time " << procMeshes()[0].facesInstance() << nl
                << "Reconstruction must start from this earlier time"
                << exit(FatalError);
        }

        reconstruct();
    }

    return !load;
}


void Foam::domainDecomposition::postReadReconstruct()
{
    forAll(procMeshes_, proci)
    {
        procMeshes_[proci].stitcher().connect(false, false, true);
    }
}


void Foam::domainDecomposition::unconformReadReconstruct()
{
    if (!procsConformal())
    {
        procFaceAddressingBf_.clear();

        completeMesh_->conform();

        unconform();
    }
}


void Foam::domainDecomposition::writeReadReconstruct
(
    const bool reconstructed,
    const bool doSets
)
{
    writeComplete(doSets);

    if (reconstructed)
    {
        writeCompletePoints(procMeshes()[0].facesInstance());
    }
}


Foam::fvMesh::readUpdateState Foam::domainDecomposition::readUpdateDecompose()
{
    const fvMesh::readUpdateState stat = readUpdate();

    // Topology changes
    {
        const label facesCompare =
            compareInstances
            (
                completeMesh().facesInstance(),
                procMeshes_[0].facesInstance()
            );

        // If the complete mesh has newer topology then we need to decompose
        if (facesCompare == -1)
        {
            decompose();
        }

        // If there has been matching topology change then reload the addressing
        if (facesCompare == 0 && stat >= fvMesh::TOPO_CHANGE)
        {
            readAddressing();
        }

        // The processor meshes should not have newer topology when decomposing
        if (facesCompare == +1)
        {
            FatalErrorInFunction
                << "Cannot decompose at time "
                << procMeshes_[0].facesInstance()
                << " because the processor mesh topology has evolved further"
                << " than the complete mesh topology." << exit(FatalError);
        }
    }

    // Geometry changes
    {
        const label pointsCompare =
            compareInstances
            (
                completeMesh().pointsInstance(),
                procMeshes_[0].pointsInstance()
            );

        // If the complete mesh has newer geometry then we need to decompose
        // the points
        if (pointsCompare == -1)
        {
            decomposePoints();
        }

        // The processor meshes should not have newer geometry when decomposing
        if (pointsCompare == +1)
        {
            FatalErrorInFunction
                << "Cannot decompose at time "
                << procMeshes_[0].pointsInstance()
                << " because the processor mesh geometry has evolved further"
                << " than the complete mesh geometry." << exit(FatalError);
        }
    }

    return stat;
}


void Foam::domainDecomposition::postReadUpdateDecompose
(
    const fvMesh::readUpdateState stat
)
{
    if (completeMesh_->stitcher().stitches() && stat != fvMesh::UNCHANGED)
    {
        procFaceAddressingBf_.clear();

        completeMesh_->stitcher().connect(false, false, true);
    }
}


void Foam::domainDecomposition::unconformReadUpdateDecompose
(
    const fvMesh::readUpdateState stat
)
{
    if (completeMesh_->stitcher().stitches() && stat != fvMesh::UNCHANGED)
    {
        forAll(procMeshes_, proci)
        {
            procMeshes_[proci].conform();
        }

        unconform();
    }
}


Foam::fvMesh::readUpdateState Foam::domainDecomposition::readUpdateReconstruct()
{
    const fvMesh::readUpdateState stat = readUpdate();

    // Topology changes
    {
        const label facesCompare =
            compareInstances
            (
                completeMesh().facesInstance(),
                procMeshes_[0].facesInstance()
            );

        // The complete mesh should not have newer topology when reconstructing
        if (facesCompare == -1)
        {
            FatalErrorInFunction
                << "Cannot reconstruct at time "
                << completeMesh().facesInstance()
                << " because the complete mesh topology has evolved further"
                << " than the processor mesh topology." << exit(FatalError);
        }

        // If there has been matching topology change then reload the addressing
        if (facesCompare == 0 && stat >= fvMesh::TOPO_CHANGE)
        {
            readAddressing();
        }

        // If the processor meshes have newer topology then we need to
        // reconstruct
        if (facesCompare == +1)
        {
            reconstruct();
        }
    }

    // Geometry changes
    {
        const label pointsCompare =
            compareInstances
            (
                completeMesh().pointsInstance(),
                procMeshes_[0].pointsInstance()
            );

        // The complete mesh should not have newer geometry when reconstructing
        if (pointsCompare == -1)
        {
            FatalErrorInFunction
                << "Cannot reconstruct at time "
                << completeMesh().pointsInstance()
                << " because the complete mesh geometry has evolved further"
                << " than the processor mesh geometry." << exit(FatalError);
        }

        // If the processor meshes have newer geometry then we need to
        // reconstruct the points
        if (pointsCompare == +1)
        {
            reconstructPoints();
        }
    }

    return stat;
}


void Foam::domainDecomposition::postReadUpdateReconstruct
(
    const fvMesh::readUpdateState stat
)
{
    if (completeMesh_->stitcher().stitches() && stat != fvMesh::UNCHANGED)
    {
        procFaceAddressingBf_.clear();

        forAll(procMeshes_, proci)
        {
            procMeshes_[proci].stitcher().connect(false, false, true);
        }
    }
}


void Foam::domainDecomposition::unconformReadUpdateReconstruct
(
    const fvMesh::readUpdateState stat
)
{
    if (completeMesh_->stitcher().stitches() && stat != fvMesh::UNCHANGED)
    {
        completeMesh_->conform();

        unconform();
    }
}


const Foam::PtrList<Foam::surfaceLabelField::Boundary>&
Foam::domainDecomposition::procFaceAddressingBf() const
{
    validateComplete();
    validateProcs();

    if (procFaceAddressingBf_.empty())
    {
        // Get any non-conformal proc-face addressing
        List<List<DynamicList<label>>> nonConformalProcFaceAddressingBf =
            this->nonConformalProcFaceAddressingBf();

        // Build finite volume face addressing boundary fields
        procFaceAddressingBf_.resize(nProcs());
        forAll(procMeshes_, proci)
        {
            const fvMesh& procMesh = procMeshes_[proci];

            procFaceAddressingBf_.set
            (
                proci,
                new surfaceLabelField::Boundary
                (
                    procMesh.boundary(),
                    surfaceLabelField::null(),
                    calculatedFvsPatchLabelField::typeName
                )
            );

            forAll(procMesh.boundary(), procPatchi)
            {
                const fvPatch& fvp = procMesh.boundary()[procPatchi];

                if (isA<nonConformalFvPatch>(fvp))
                {
                    procFaceAddressingBf_[proci][procPatchi] =
                        nonConformalProcFaceAddressingBf[proci][procPatchi];
                }
                else if (isA<processorCyclicFvPatch>(fvp))
                {
                    const label completePatchi =
                        refCast<const processorCyclicFvPatch>(fvp)
                       .referPatchIndex();

                    procFaceAddressingBf_[proci][procPatchi] =
                        mag(fvp.patchSlice(procFaceAddressing_[proci]))
                      - completeMesh().boundaryMesh()[completePatchi].start();
                }
                else if (isA<processorFvPatch>(fvp))
                {
                    procFaceAddressingBf_[proci][procPatchi] =
                        fvp.patchSlice(procFaceAddressing_[proci]);
                }
                else
                {
                    procFaceAddressingBf_[proci][procPatchi] =
                        mag(fvp.patchSlice(procFaceAddressing_[proci]))
                      - completeMesh().boundaryMesh()[procPatchi].start();
                }
            }
        }
    }

    return procFaceAddressingBf_;
}


void Foam::domainDecomposition::writeComplete(const bool doSets) const
{
    const bool topologyWrite =
        static_cast<const faceCompactIOList&>(completeMesh().faces())
       .writeOpt() == IOobject::AUTO_WRITE;

    // Set the precision of the points data to be min 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    // Write the complete mesh
    completeMesh().write();

    // Everything written below is topological data, so quit here if not
    // writing topology
    if (!topologyWrite) return;

    // Read, reconstruct and write sets
    if (doSets)
    {
        HashPtrTable<cellSet> cellSets;
        HashPtrTable<faceSet> faceSets;
        HashPtrTable<pointSet> pointSets;

        for (label proci = 0; proci < nProcs(); proci++)
        {
            const fvMesh& procMesh = procMeshes_[proci];

            const labelList& cellMap = procCellAddressing()[proci];
            const labelList& faceMap = procFaceAddressing()[proci];
            const labelList& pointMap = procPointAddressing()[proci];

            // Scan the contents of the sets directory
            IOobjectList setObjects
            (
                procMesh,
                runTimes_.procTimes()[proci].name(),
                polyMesh::meshSubDir/"sets"
            );
            IOobjectList cellSetObjects
            (
                setObjects.lookupClass(cellSet::typeName)
            );
            IOobjectList faceSetObjects
            (
                setObjects.lookupClass(faceSet::typeName)
            );
            IOobjectList pointSetObjects
            (
                setObjects.lookupClass(pointSet::typeName)
            );

            // Read and reconstruct the sets
            forAllConstIter(IOobjectList, cellSetObjects, iter)
            {
                const cellSet procSet(*iter());

                if (!cellSets.found(iter.key()))
                {
                    cellSets.insert
                    (
                        iter.key(),
                        new cellSet
                        (
                            completeMesh(),
                            iter.key(),
                            procSet.size()
                        )
                    );
                }

                cellSet& cSet = *cellSets[iter.key()];

                cSet.instance() = runTimes_.completeTime().name();

                forAllConstIter(cellSet, procSet, iter)
                {
                    cSet.insert(cellMap[iter.key()]);
                }
            }
            forAllConstIter(IOobjectList, faceSetObjects, iter)
            {
                const faceSet procSet(*iter());

                if (!faceSets.found(iter.key()))
                {
                    faceSets.insert
                    (
                        iter.key(),
                        new faceSet
                        (
                            completeMesh(),
                            iter.key(),
                            procSet.size()
                        )
                    );
                }

                faceSet& cSet = *faceSets[iter.key()];

                cSet.instance() = runTimes_.completeTime().name();

                forAllConstIter(faceSet, procSet, iter)
                {
                    cSet.insert(mag(faceMap[iter.key()]) - 1);
                }
            }
            forAllConstIter(IOobjectList, pointSetObjects, iter)
            {
                const pointSet procSet(*iter());

                if (!pointSets.found(iter.key()))
                {
                    pointSets.insert
                    (
                        iter.key(),
                        new pointSet
                        (
                            completeMesh(),
                            iter.key(),
                            procSet.size()
                        )
                    );
                }

                pointSet& cSet = *pointSets[iter.key()];

                cSet.instance() = runTimes_.completeTime().name();

                forAllConstIter(pointSet, procSet, iter)
                {
                    cSet.insert(pointMap[iter.key()]);
                }
            }
        }

        // Write the sets
        forAllConstIter(HashPtrTable<cellSet>, cellSets, iter)
        {
            iter()->write();
        }
        forAllConstIter(HashPtrTable<faceSet>, faceSets, iter)
        {
            iter()->write();
        }
        forAllConstIter(HashPtrTable<pointSet>, pointSets, iter)
        {
            iter()->write();
        }
    }

    // Read, decompose, and write refinement data (if any)
    UPtrList<const labelList> cellMaps(nProcs());
    UPtrList<const labelList> pointMaps(nProcs());
    PtrList<const hexRef8Data> refinementDatas(nProcs());
    for (label proci = 0; proci < nProcs(); proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        cellMaps.set(proci, &procCellAddressing_[proci]);
        pointMaps.set(proci, &procPointAddressing_[proci]);
        refinementDatas.set
        (
            proci,
            new hexRef8Data
            (
                IOobject
                (
                    "dummy",
                    completeMesh().facesInstance(),
                    polyMesh::meshSubDir,
                    procMesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
    }
    hexRef8Data
    (
        IOobject
        (
            "dummy",
            completeMesh().facesInstance(),
            polyMesh::meshSubDir,
            completeMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        cellMaps,
        pointMaps,
        refinementDatas
    ).write();

    // Write decomposition addressing
    writeAddressing();
}


void Foam::domainDecomposition::writeProcs(const bool doSets) const
{
    const bool topologyWrite =
        static_cast<const faceCompactIOList&>(procMeshes()[0].faces())
       .writeOpt() == IOobject::AUTO_WRITE;

    // Write out the meshes
    for (label proci = 0; proci < nProcs(); proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        // Set the precision of the points data to be min 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        // Write the processor mesh
        procMesh.write();
    }

    // Everything written below is topological data, so quit here if not
    // writing topology
    if (!topologyWrite) return;

    // Read, decompose, and write any sets
    if (doSets)
    {
        // Scan the contents of the sets directory
        IOobjectList setObjects
        (
            completeMesh(),
            completeMesh().facesInstance(),
            polyMesh::meshSubDir/"sets"
        );
        IOobjectList cellSetObjects
        (
            setObjects.lookupClass(cellSet::typeName)
        );
        IOobjectList faceSetObjects
        (
            setObjects.lookupClass(faceSet::typeName)
        );
        IOobjectList pointSetObjects
        (
            setObjects.lookupClass(pointSet::typeName)
        );

        // Read the sets
        PtrList<const cellSet> cellSets;
        forAllConstIter(IOobjectList, cellSetObjects, iter)
        {
            cellSets.append(new cellSet(*iter()));
        }
        PtrList<const faceSet> faceSets;
        forAllConstIter(IOobjectList, faceSetObjects, iter)
        {
            faceSets.append(new faceSet(*iter()));
        }
        PtrList<const pointSet> pointSets;
        forAllConstIter(IOobjectList, pointSetObjects, iter)
        {
            pointSets.append(new pointSet(*iter()));
        }

        // Decompose and write sets into the processor mesh directories
        for (label proci = 0; proci < nProcs(); proci++)
        {
            const fvMesh& procMesh = procMeshes_[proci];

            forAll(cellSets, i)
            {
                const cellSet& cs = cellSets[i];
                cellSet set(procMesh, cs.name(), cs.size()/nProcs());
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
                faceSet set(procMesh, cs.name(), cs.size()/nProcs());
                forAll(procFaceAddressing_[proci], i)
                {
                    if (cs.found(mag(procFaceAddressing_[proci][i]) - 1))
                    {
                        set.insert(i);
                    }
                }
                set.write();
            }
            forAll(pointSets, i)
            {
                const pointSet& cs = pointSets[i];
                pointSet set(procMesh, cs.name(), cs.size()/nProcs());
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
    }

    // Read, decompose, and write refinement data (if any)
    const hexRef8Data refinementData
    (
        IOobject
        (
            "dummy",
            completeMesh().facesInstance(),
            polyMesh::meshSubDir,
            completeMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );
    for (label proci = 0; proci < nProcs(); proci++)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        hexRef8Data
        (
            IOobject
            (
                "dummy",
                completeMesh().facesInstance(),
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

    // Write decomposition addressing
    writeAddressing();
}


// ************************************************************************* //
