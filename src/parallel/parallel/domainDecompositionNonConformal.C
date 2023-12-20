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
#include "cyclicFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "nonConformalCyclicFvPatch.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "nonConformalErrorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void checkNonConformalCoupledPatchOrdering
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
}


void checkNonConformalErrorPatchOrdering
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
}


}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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


Foam::List<Foam::List<Foam::DynamicList<Foam::label>>>
Foam::domainDecomposition::nonConformalProcFaceAddressingBf() const
{
    validateComplete();
    validateProcs();

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

                refPatchProcPatchTable[pcFvp.referPatchIndex()].insert
                (
                    labelPair(proci, pcFvp.neighbProcNo()),
                    procPatchi
                );
            }
        }
    }

    // Build non-conformal finite volume face addressing for each processor
    List<List<DynamicList<label>>> result(nProcs());
    forAll(result, proci)
    {
        result[proci].resize
        (
            procMeshes_[proci].boundary().size()
        );
    }

    if (completeConformal() && procsConformal())
    {
        // Nothing to do
    }
    else if (!completeConformal())
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

            const label nccNbrPatchi = nccFvp.nbrPatchIndex();

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

                result[proci][procNccPatchi].append(nccPatchFacei + 1);
                result[nbrProci][nbrProcNccPatchi].append(nccPatchFacei + 1);
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

                result[proci][ncePatchi].append(ncePatchFacei + 1);
            }
        }
    }
    else // if (!procsConformal())
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

            const label nccNbrPatchi = nccFvp.nbrPatchIndex();

            label nccPatchFacei = 0;
            labelPairTable procNccPatchFaceis;
            forAllConstIter
            (
                labelPairTable,
                refPatchProcPatchTable[nccPatchi],
                iter
            )
            {
                procNccPatchFaceis.insert(iter.key(), 0);
            }

            while (true)
            {
                labelPair procNbrProc(labelMax, labelMax);
                labelPair faceNbrFace(labelMax, labelMax);

                forAllConstIter(labelPairTable, procNccPatchFaceis, iter)
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

                    result[proci][procNccPatchi]
                        .append(nccPatchFacei + 1);
                    result[nbrProci][nbrProcNccPatchi]
                        .append(nccPatchFacei + 1);

                    nccPatchFacei ++;
                    procNccPatchFaceis[procNbrProc] ++;
                }
            }
        }

        // Error patches
        forAll(completeMesh().boundary(), ncePatchi)
        {
            const fvPatch& fvp = completeMesh().boundary()[ncePatchi];

            if (!isA<nonConformalErrorFvPatch>(fvp)) continue;

            label ncePatchFacei = 0;
            labelList procNcePatchFaceis(nProcs(), 0);

            while (true)
            {
                label facei = labelMax, proci = labelMax;

                forAll(procNcePatchFaceis, procStari)
                {
                    const label size =
                        procMeshes_[procStari]
                       .polyFacesBf()[ncePatchi]
                       .size();

                    if (procNcePatchFaceis[procStari] >= size) continue;

                    const label procFacei =
                        procMeshes_[procStari].polyFacesBf()
                        [ncePatchi][procNcePatchFaceis[procStari]];

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
                    result[proci][ncePatchi].append(ncePatchFacei + 1);

                    ncePatchFacei ++;
                    procNcePatchFaceis[proci] ++;
                }
            }
        }
    }

    return result;
}


void Foam::domainDecomposition::unconformComplete()
{
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
                   .referPatchIndex()
                  : procNccPatchi;

                const label size =
                    max
                    (
                        max(mag(faceAddressingBf[procNccPatchi])),
                        polyFacesBf[completeNccPatchi].size()
                    );

                polyFacesBf[completeNccPatchi].resize(size, -1);
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

                // Set dummy data for the face geometry. This should not be
                // used during decomposition.
                Sf.boundaryFieldRef()[completeNccPatchi].resize(size, Zero);
                Cf.boundaryFieldRef()[completeNccPatchi].resize(size, Zero);
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
                   .nbrPatchIndex();

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


void Foam::domainDecomposition::unconformProcs()
{
    // Construct the reverse of proc-face-face addressing. -1 indicates a face
    // that is on a (conformal) processor boundary and hence has multiple
    // associated proc-face indices.
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
                   .referPatchIndex()
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

                // Set dummy data for the face geometry. This should not be
                // used during decomposition.
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
                            refCast<const cyclicFvPatch>(fvp).nbrPatchIndex();
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
                                 && nbrPcFvp.referPatchIndex()
                                 == pcFvp.referPatch().nbrPatchIndex()
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


void Foam::domainDecomposition::unconform()
{
    if (completeConformal() && procsConformal())
    {
        // Nothing to do
    }
    else if (!completeConformal() && procsConformal())
    {
        unconformProcs();
    }
    else if (completeConformal() && !procsConformal())
    {
        unconformComplete();
    }
    else // if (!completeConformal() && !procsConformal())
    {
        // Assume everything is consistent and do nothing
    }
}


// ************************************************************************* //
