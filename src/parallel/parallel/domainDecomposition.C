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
#include "decompositionMethod.H"
#include "IOobjectList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "hexRef8Data.H"
#include "cyclicFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "nonConformalCyclicFvPatch.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "nonConformalErrorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(domainDecomposition, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::domainDecomposition::decomposePoints()
{
    for (label proci = 0; proci < nProcs(); proci++)
    {
        fvMesh& procMesh = procMeshes_[proci];

        const label pointsCompare =
            compareInstances
            (
                completeMesh().pointsInstance(),
                procMeshes_[proci].pointsInstance()
            );

        if (pointsCompare == -1)
        {
            procMesh.setPoints
            (
                pointField
                (
                    completeMesh().points(),
                    procPointAddressing_[proci]
                )
            );
        }
    }
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


bool Foam::domainDecomposition::completeConformal() const
{
    return completeMesh().conformal();
}


bool Foam::domainDecomposition::procsConformal() const
{
    forAll(procMeshes_, proci)
    {
        if (!procMeshes_[proci].conformal())
        {
            return false;
        }
    }

    return true;
}


void Foam::domainDecomposition::unconform()
{
    // ...
    auto checkNonConformalCoupledPatchOrdering = []
    (
        const labelPair& procs,
        const fvPatch& fvp,
        const fvPatch& nbrFvp,
        const labelList& polyFacesPf,
        const labelList& nbrPolyFacesPf
    )
    {
        if (fvp.size() != nbrFvp.size())
        {
            FatalErrorInFunction
                << "Coupled patches " << fvp.name() << " and "
                << nbrFvp.name() << " are not the same size"
                << exit(FatalError);
        }

        if (fvp.size())
        {
            for (label i = 1; i < fvp.size(); ++ i)
            {
                if
                (
                    polyFacesPf[i - 1] > polyFacesPf[i]
                  ? true
                  : polyFacesPf[i - 1] == polyFacesPf[i]
                  ? nbrPolyFacesPf[i - 1] >= nbrPolyFacesPf[i]
                  : false
                )
                {
                    FatalErrorInFunction
                        << "Coupled patches " << fvp.name()
                        << " and " << nbrFvp.name()
                        << " are not in order";

                    if (procs[0] == procs[1])
                    {
                        FatalErrorInFunction
                            << " on process #" << procs[0];
                    }

                    FatalErrorInFunction
                        << exit(FatalError);
                }
            }
        }

        InfoInFunction
            << "Coupled patches " << fvp.name()
            << " and " << nbrFvp.name()
            << " are in order" << endl;
    };

    // ...
    auto checkNonConformalErrorPatchOrdering = []
    (
        const label& proci,
        const fvPatch& fvp,
        const labelList& polyFacesPf
    )
    {
        if (fvp.size())
        {
            for (label i = 1; i < fvp.size(); ++ i)
            {
                if (polyFacesPf[i - 1] > polyFacesPf[i])
                {
                    FatalErrorInFunction
                        << "Error patch " << fvp.name()
                        << " is not in order";

                    if (proci > 0)
                    {
                        FatalErrorInFunction
                            << " on process #" << proci;
                    }

                    FatalErrorInFunction
                        << exit(FatalError);
                }
            }
        }

        InfoInFunction
            << "Error patch " << fvp.name()
            << " is in order" << endl;
    };

    if (completeConformal() && procsConformal())
    {
        // Nothing to do
    }
    else if (!completeConformal() && procsConformal())
    {
        // Construct the reverse of proc-face-face addressing
        labelList faceProcFace(completeMesh().nFaces(), -labelMax);
        forAll(procMeshes_, proci)
        {
            forAll(procFaceAddressing()[proci], procFacei)
            {
                const label facei =
                    mag(procFaceAddressing()[proci][procFacei]) - 1;

                faceProcFace[facei] =
                    faceProcFace[facei] == -labelMax ? procFacei : -1;
            }
        }

        // Map the polyFacesBf from the complete to the processor meshes and
        // unconform the processor meshes. Set dummy data for the face
        // geometry. This should not be used during decomposition.
        forAll(procMeshes_, proci)
        {
            fvMesh& procMesh = procMeshes_[proci];

            const surfaceLabelField::Boundary& faceAddressingBf =
                procFaceAddressingBf()[proci];

            surfaceLabelField::Boundary polyFacesBf
            (
                surfaceLabelField::null(),
                procMesh.polyFacesBf()
            );
            surfaceVectorField Sf(procMesh.Sf().cloneUnSliced());
            surfaceVectorField Cf(procMesh.Cf().cloneUnSliced());

            forAll(procMesh.boundary(), procNccPatchi)
            {
                const fvPatch& fvp = procMesh.boundary()[procNccPatchi];

                if (isA<nonConformalFvPatch>(fvp))
                {
                    const label completeNccPatchi =
                        isA<processorCyclicFvPatch>(fvp)
                      ? refCast<const processorCyclicFvPatch>(fvp)
                       .referPatchID()
                      : procNccPatchi;

                    polyFacesBf[procNccPatchi] =
                        labelField
                        (
                            faceProcFace,
                            labelField
                            (
                                completeMesh().polyFacesBf()[completeNccPatchi],
                                mag(faceAddressingBf[procNccPatchi]) - 1
                            )
                        );

                    const label size = polyFacesBf[procNccPatchi].size();

                    Sf.boundaryFieldRef()[procNccPatchi].resize(size, Zero);
                    Cf.boundaryFieldRef()[procNccPatchi].resize(size, Zero);
                }
            }

            procMesh.unconform
            (
                polyFacesBf,
                Sf,
                Cf,
                NullObjectRef<surfaceScalarField>(),
                false
            );
        }

        // Check ordering
        if (debug)
        {
            forAll(procMeshes_, proci)
            {
                forAll(procMeshes_[proci].boundary(), patchi)
                {
                    const fvPatch& fvp =
                        procMeshes_[proci].boundary()[patchi];

                    // Coupled patches
                    if
                    (
                        isA<nonConformalCoupledFvPatch>(fvp)
                     && refCast<const nonConformalCoupledFvPatch>(fvp).owner()
                    )
                    {
                        const label nccPatchi = patchi;

                        label nbrProci = -1, nbrNccPatchi = -1;
                        if (isA<cyclicFvPatch>(fvp))
                        {
                            nbrProci = proci;
                            nbrNccPatchi =
                                refCast<const cyclicFvPatch>(fvp).nbrPatchID();
                        }
                        else if (isA<processorCyclicFvPatch>(fvp))
                        {
                            typedef processorCyclicFvPatch PcFvp;

                            const PcFvp& pcFvp = refCast<const PcFvp>(fvp);

                            nbrProci = pcFvp.neighbProcNo();

                            const fvBoundaryMesh& nbrFvPatches =
                                procMeshes_[nbrProci].boundary();

                            forAll(nbrFvPatches, nbrNccPatchj)
                            {
                                const fvPatch& nbrFvp =
                                    nbrFvPatches[nbrNccPatchj];

                                if (isA<PcFvp>(nbrFvp))
                                {
                                    const PcFvp& nbrPcFvp =
                                        refCast<const PcFvp>(nbrFvp);

                                    if
                                    (
                                        nbrPcFvp.neighbProcNo()
                                     == proci
                                     && nbrPcFvp.referPatchID()
                                     == pcFvp.referPatch().nbrPatchID()
                                    )
                                    {
                                        nbrNccPatchi = nbrNccPatchj;
                                        break;
                                    }
                                }
                            }

                            if (nbrNccPatchi == -1)
                            {
                                FatalErrorInFunction
                                    << "Opposite processor patch not found for "
                                    << "patch " << fvp.name() << " on proc #"
                                    << proci << exit(FatalError);
                            }
                        }
                        else
                        {
                            FatalErrorInFunction
                                << "Non-conformal-coupled type not recognised "
                                << "for patch " << fvp.name() << " on proc #"
                                << proci << exit(FatalError);
                        }

                        checkNonConformalCoupledPatchOrdering
                        (
                            {proci, nbrProci},
                            procMeshes_[proci].boundary()[nccPatchi],
                            procMeshes_[nbrProci].boundary()[nbrNccPatchi],
                            procMeshes_[proci].polyFacesBf()[nccPatchi],
                            procMeshes_[nbrProci].polyFacesBf()[nbrNccPatchi]
                        );
                    }

                    // Error patches
                    if (isA<nonConformalErrorFvPatch>(fvp))
                    {
                        const label ncePatchi = patchi;

                        checkNonConformalErrorPatchOrdering
                        (
                            proci,
                            procMeshes_[proci].boundary()[ncePatchi],
                            procMeshes_[proci].polyFacesBf()[ncePatchi]
                        );
                    }
                }
            }
        }
    }
    else if (completeConformal() && !procsConformal())
    {
        // Map the polyFacesBf from the processor to the complete meshes and
        // unconform the complete mesh. Set dummy data for the face
        // geometry. This should not be used during decomposition.
        surfaceLabelField::Boundary polyFacesBf
        (
            surfaceLabelField::null(),
            completeMesh().polyFacesBf()
        );
        surfaceVectorField Sf(completeMesh().Sf().cloneUnSliced());
        surfaceVectorField Cf(completeMesh().Cf().cloneUnSliced());

        forAll(procMeshes_, proci)
        {
            const fvMesh& procMesh = procMeshes_[proci];

            const surfaceLabelField::Boundary& faceAddressingBf =
                procFaceAddressingBf()[proci];

            forAll(procMesh.boundary(), procNccPatchi)
            {
                const fvPatch& fvp = procMesh.boundary()[procNccPatchi];

                if (isA<nonConformalFvPatch>(fvp))
                {
                    const label completeNccPatchi =
                        isA<processorCyclicFvPatch>(fvp)
                      ? refCast<const processorCyclicFvPatch>(fvp)
                       .referPatchID()
                      : procNccPatchi;

                    const label size =
                        max
                        (
                            max(mag(faceAddressingBf[procNccPatchi])),
                            polyFacesBf[completeNccPatchi].size()
                        );

                    polyFacesBf[completeNccPatchi].resize(size, -1);
                    Sf.boundaryFieldRef()[completeNccPatchi].resize(size, Zero);
                    Cf.boundaryFieldRef()[completeNccPatchi].resize(size, Zero);

                    polyFacesBf[completeNccPatchi].labelField::rmap
                    (
                        mag
                        (
                            labelField
                            (
                                procFaceAddressing_[proci],
                                procMesh.polyFacesBf()[procNccPatchi]
                            )
                        ) - 1,
                        mag(faceAddressingBf[procNccPatchi]) - 1
                    );
                }
            }
        }

        completeMesh_->unconform
        (
            polyFacesBf,
            Sf,
            Cf,
            NullObjectRef<surfaceScalarField>(),
            false
        );

        // Check ordering
        if (debug)
        {
            forAll(completeMesh().boundary(), patchi)
            {
                const fvPatch& fvp = completeMesh().boundary()[patchi];

                // Coupled patches
                if
                (
                    isA<nonConformalCyclicFvPatch>(fvp)
                 && refCast<const nonConformalCoupledFvPatch>(fvp).owner()
                )
                {
                    const label nccPatchi = patchi;

                    const label nbrNccPatchi =
                        refCast<const nonConformalCyclicFvPatch>(fvp)
                       .nbrPatchID();

                    checkNonConformalCoupledPatchOrdering
                    (
                        {-labelMax, labelMax},
                        completeMesh().boundary()[nccPatchi],
                        completeMesh().boundary()[nbrNccPatchi],
                        completeMesh().polyFacesBf()[nccPatchi],
                        completeMesh().polyFacesBf()[nbrNccPatchi]
                    );
                }

                // Error patches
                if (isA<nonConformalErrorFvPatch>(fvp))
                {
                    const label ncePatchi = patchi;

                    checkNonConformalErrorPatchOrdering
                    (
                        -labelMax,
                        completeMesh().boundary()[ncePatchi],
                        completeMesh().polyFacesBf()[ncePatchi]
                    );
                }
            }
        }
    }
    else // if (!completeConformal() && !procsConformal())
    {
        // Assume everything is consistent and do nothing
    }
}


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


void Foam::domainDecomposition::readComplete(const bool stitch)
{
    completeMesh_.reset
    (
        new fvMesh
        (
            IOobject
            (
                regionName_,
                runTimes_.completeTime().name(),
                runTimes_.completeTime(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            false,
            stitch
          ? fvMesh::stitchType::nonGeometric
          : fvMesh::stitchType::none
        )
    );
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
                    runTimes_.procTimes()[proci],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                false
            )
        );
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
}


Foam::fvMesh::readUpdateState Foam::domainDecomposition::readUpdate()
{
    validateComplete();
    validateProcs();

    // Do read-update on all meshes
    fvMesh::readUpdateState stat =
        completeMesh_->readUpdate(fvMesh::stitchType::nonGeometric);
    forAll(runTimes_.procTimes(), proci)
    {
        fvMesh::readUpdateState procStat =
            procMeshes_[proci].readUpdate(fvMesh::stitchType::nonGeometric);
        if (procStat > stat)
        {
            stat = procStat;
        }
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
    const word& regionName
)
:
    runTimes_(runTimes),
    regionName_(regionName),
    completeMesh_(nullptr),
    procMeshes_(nProcs()),
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

bool Foam::domainDecomposition::readDecompose(const bool doSets)
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

    if (!completeConformal())
    {
        procFaceAddressingBf_.clear();
        forAll(procMeshes_, proci) procMeshes_[proci].conform();
        unconform();
    }

    writeProcs(doSets);

    if (!load)
    {
        writeProcPoints(completeMesh().facesInstance());
    }

    return !load;
}


bool Foam::domainDecomposition::readReconstruct(const bool doSets)
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

        readComplete(completePointsIo.headerOk());

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

    if (!procsConformal())
    {
        procFaceAddressingBf_.clear();
        completeMesh_->conform();
        unconform();
    }

    writeComplete(doSets);

    if (!load)
    {
        writeCompletePoints(procMeshes()[0].facesInstance());
    }

    return !load;
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
            procFaceAddressingBf_.clear();
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

    // Non-conformal changes
    {
        // If the mesh has changed in any way, and the complete mesh is
        // non-conformal, then we need to re-unconform the processor meshes
        if (stat != fvMesh::UNCHANGED && !completeConformal())
        {
            procFaceAddressingBf_.clear();
            forAll(procMeshes_, proci) procMeshes_[proci].conform();
            unconform();
        }
    }

    return stat;
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
            procFaceAddressingBf_.clear();
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

    // Non-conformal changes
    {
        // If the mesh has changed in any way, and the processor meshes are
        // non-conformal, then we need to re-unconform the complete mesh
        if (stat != fvMesh::UNCHANGED && !procsConformal())
        {
            procFaceAddressingBf_.clear();
            completeMesh_->conform();
            unconform();
        }
    }

    return stat;
}


const Foam::PtrList<Foam::surfaceLabelField::Boundary>&
Foam::domainDecomposition::procFaceAddressingBf() const
{
    validateComplete();
    validateProcs();

    if (procFaceAddressingBf_.empty())
    {
        // Map from reference patch and processors to the interface patch
        typedef HashTable<label, labelPair, Hash<labelPair>> labelPairTable;
        List<labelPairTable> refPatchProcPatchTable
        (
            completeMesh().boundary().size()
        );
        forAll(procMeshes_, proci)
        {
            const fvMesh& procMesh = procMeshes_[proci];

            forAll(procMesh.boundary(), procPatchi)
            {
                const fvPatch& fvp = procMesh.boundary()[procPatchi];

                if (isA<cyclicFvPatch>(fvp))
                {
                    refPatchProcPatchTable[procPatchi].insert
                    (
                        labelPair(proci, proci),
                        procPatchi
                    );
                }
                else if (isA<processorCyclicFvPatch>(fvp))
                {
                    const processorCyclicFvPatch& pcFvp =
                        refCast<const processorCyclicFvPatch>(fvp);

                    refPatchProcPatchTable[pcFvp.referPatchID()].insert
                    (
                        labelPair(proci, pcFvp.neighbProcNo()),
                        procPatchi
                    );
                }
            }
        }

        // Build non-conformal finite volume face addressing for each processor
        List<List<DynamicList<label>>>
            nonConformalProcFaceAddressingBf(nProcs());
        forAll(nonConformalProcFaceAddressingBf, proci)
        {
            nonConformalProcFaceAddressingBf[proci].resize
            (
                procMeshes_[proci].boundary().size()
            );
        }
        if (completeMesh().conformal() && procMeshes_[0].conformal())
        {
            // Nothing to do
        }
        else if (!completeMesh().conformal())
        {
            // Decompose non-conformal addressing
            const surfaceLabelField::Boundary& polyFacesBf =
                completeMesh().polyFacesBf();

            // Cyclic patches
            forAll(completeMesh().boundary(), nccPatchi)
            {
                const fvPatch& fvp = completeMesh().boundary()[nccPatchi];

                if (!isA<nonConformalCyclicFvPatch>(fvp)) continue;

                const nonConformalCyclicFvPatch& nccFvp =
                    refCast<const nonConformalCyclicFvPatch>(fvp);

                if (!nccFvp.owner()) continue;

                const label nccNbrPatchi = nccFvp.nbrPatchID();

                forAll(polyFacesBf[nccPatchi], nccPatchFacei)
                {
                    const label facei = polyFacesBf[nccPatchi][nccPatchFacei];
                    const label celli = completeMesh().faceOwner()[facei];
                    const label proci = cellProc_[celli];

                    const label nbrFacei =
                        polyFacesBf[nccNbrPatchi][nccPatchFacei];
                    const label nbrCelli =
                        completeMesh().faceOwner()[nbrFacei];
                    const label nbrProci = cellProc_[nbrCelli];

                    const label procNccPatchi =
                        refPatchProcPatchTable
                        [nccPatchi][labelPair(proci, nbrProci)];
                    const label nbrProcNccPatchi =
                        refPatchProcPatchTable
                        [nccNbrPatchi][labelPair(nbrProci, proci)];

                    nonConformalProcFaceAddressingBf[proci][procNccPatchi]
                       .append(nccPatchFacei + 1);
                    nonConformalProcFaceAddressingBf[nbrProci][nbrProcNccPatchi]
                       .append(nccPatchFacei + 1);
                }
            }

            // Error patches
            forAll(completeMesh().boundary(), ncePatchi)
            {
                const fvPatch& fvp = completeMesh().boundary()[ncePatchi];

                if (!isA<nonConformalErrorFvPatch>(fvp)) continue;

                forAll(polyFacesBf[ncePatchi], ncePatchFacei)
                {
                    const label facei = polyFacesBf[ncePatchi][ncePatchFacei];
                    const label celli = completeMesh().faceOwner()[facei];
                    const label proci = cellProc_[celli];

                    nonConformalProcFaceAddressingBf[proci][ncePatchi]
                       .append(ncePatchFacei + 1);
                }
            }
        }
        else // if (!procMeshes_[0].conformal())
        {
            // Reconstruct non-conformal addressing

            // Cyclic patches
            forAll(completeMesh().boundary(), nccPatchi)
            {
                const fvPatch& fvp = completeMesh().boundary()[nccPatchi];

                if (!isA<nonConformalCyclicFvPatch>(fvp)) continue;

                const nonConformalCyclicFvPatch& nccFvp =
                    refCast<const nonConformalCyclicFvPatch>(fvp);

                if (!nccFvp.owner()) continue;

                const label nccNbrPatchi = nccFvp.nbrPatchID();

                label count = 0;
                labelPairTable procCounts;
                forAllConstIter
                (
                    labelPairTable,
                    refPatchProcPatchTable[nccPatchi],
                    iter
                )
                {
                    procCounts.insert(iter.key(), 0);
                }

                while (true)
                {
                    labelPair procNbrProc(labelMax, labelMax);
                    labelPair faceNbrFace(labelMax, labelMax);

                    forAllConstIter(labelPairTable, procCounts, iter)
                    {
                        const label proci = iter.key().first();
                        const label nbrProci = iter.key().second();

                        const labelPair procNbrProcStar(proci, nbrProci);
                        const labelPair nbrProcProcStar(nbrProci, proci);

                        const label procNccPatchi =
                            refPatchProcPatchTable
                            [nccPatchi][procNbrProcStar];
                        const label nbrProcNccPatchi =
                            refPatchProcPatchTable
                            [nccNbrPatchi][nbrProcProcStar];

                        const label size =
                            procMeshes_[proci]
                           .polyFacesBf()[procNccPatchi]
                           .size();

                        if (iter() >= size) continue;

                        const label procFacei =
                            procMeshes_[proci].polyFacesBf()
                            [procNccPatchi][iter()];
                        const label nbrProcFacei =
                            procMeshes_[nbrProci].polyFacesBf()
                            [nbrProcNccPatchi][iter()];

                        const labelPair faceNbrFaceStar
                        (
                            procFaceAddressing_[proci][procFacei] - 1,
                            procFaceAddressing_[nbrProci][nbrProcFacei] - 1
                        );

                        if (faceNbrFace > faceNbrFaceStar)
                        {
                            procNbrProc = procNbrProcStar;
                            faceNbrFace = faceNbrFaceStar;
                        }
                    }

                    if (faceNbrFace == labelPair(labelMax, labelMax))
                    {
                        break;
                    }
                    else
                    {
                        const label proci = procNbrProc.first();
                        const label nbrProci = procNbrProc.second();

                        const labelPair nbrProcProc(nbrProci, proci);

                        const label procNccPatchi =
                            refPatchProcPatchTable[nccPatchi][procNbrProc];
                        const label nbrProcNccPatchi =
                            refPatchProcPatchTable[nccNbrPatchi][nbrProcProc];

                        nonConformalProcFaceAddressingBf
                            [proci][procNccPatchi].append(count + 1);
                        nonConformalProcFaceAddressingBf
                            [nbrProci][nbrProcNccPatchi].append(count + 1);

                        count ++;
                        procCounts[procNbrProc] ++;
                    }
                }
            }

            // Error patches
            forAll(completeMesh().boundary(), ncePatchi)
            {
                const fvPatch& fvp = completeMesh().boundary()[ncePatchi];

                if (!isA<nonConformalErrorFvPatch>(fvp)) continue;

                label count = 0;
                labelList procCounts(nProcs(), 0);

                while (true)
                {
                    label facei = labelMax, proci = labelMax;

                    forAll(procCounts, procStari)
                    {
                        const label size =
                            procMeshes_[procStari]
                           .polyFacesBf()[ncePatchi]
                           .size();

                        if (procCounts[procStari] >= size) continue;

                        const label procFacei =
                            procMeshes_[procStari].polyFacesBf()
                            [ncePatchi][procCounts[procStari]];

                        const label faceStari =
                            procFaceAddressing_[procStari][procFacei] - 1;

                        if (facei > faceStari)
                        {
                            facei = faceStari;
                            proci = procStari;
                        }
                    }

                    if (facei == labelMax)
                    {
                        break;
                    }
                    else
                    {
                        nonConformalProcFaceAddressingBf
                            [proci][ncePatchi].append(count + 1);

                        count ++;
                        procCounts[proci] ++;
                    }
                }
            }
        }

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
                       .referPatchID();

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

            if
            (
                (cellSets.empty() && !cellSetObjects.empty())
             || (faceSets.empty() && !faceSetObjects.empty())
             || (pointSets.empty() && !pointSetObjects.empty())
            )
            {
                Info<< "Reconstructing sets" << incrIndent << nl << endl;
            }

            // Read and reconstruct the sets
            forAllConstIter(IOobjectList, cellSetObjects, iter)
            {
                const cellSet procSet(*iter());

                if (!cellSets.found(iter.key()))
                {
                    Info<< indent << "cellSet " << iter.key() << endl;

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
                    Info<< indent << "faceSet " << iter.key() << endl;

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
                    Info<< indent << "pointSet " << iter.key() << endl;

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

        if (!cellSets.empty() || !faceSets.empty() || !pointSets.empty())
        {
            Info<< decrIndent << endl;
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
