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
#include "cyclicFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "nonConformalCyclicFvPatch.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "nonConformalMappedWallFvPatch.H"
#include "nonConformalErrorFvPatch.H"
#include "multiDomainDecomposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void checkNonConformalCoupledPatchOrdering
(
    const labelPair& procs,
    const fvPatch& fvp,
    const fvPatch& nbrFvp,
    const labelUList& polyFacesPf,
    const labelUList& nbrPolyFacesPf
)
{
    if (polyFacesPf.size() != nbrPolyFacesPf.size())
    {
        FatalErrorInFunction
            << "Coupled patches " << fvp.name() << " and "
            << nbrFvp.name() << " are not the same size"
            << exit(FatalError);
    }

    if (polyFacesPf.size())
    {
        for (label i = 1; i < polyFacesPf.size(); ++ i)
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


void checkCompleteMeshOrdering
(
    const fvMesh& completeMesh,
    const multiDomainDecomposition& regionMeshes
)
{
    forAll(completeMesh.boundary(), patchi)
    {
        const fvPatch& fvp = completeMesh.boundary()[patchi];

        // Cyclic patches
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
                completeMesh.boundary()[nccPatchi],
                completeMesh.boundary()[nbrNccPatchi],
                completeMesh.polyFacesBf()[nccPatchi],
                completeMesh.polyFacesBf()[nbrNccPatchi]
            );
        }

        // Mapped patches
        if
        (
            isA<nonConformalMappedWallFvPatch>(fvp)
         && refCast<const nonConformalMappedWallFvPatch>(fvp).owner()
        )
        {
            const nonConformalMappedWallFvPatch& ncmwFvp =
                refCast<const nonConformalMappedWallFvPatch>(fvp);

            const domainDecomposition& nbrDecomposition =
                regionMeshes[ncmwFvp.nbrRegionName()]();

            const fvMesh& nbrCompleteMesh = nbrDecomposition.completeMesh();

            const label ncmwPatchi = patchi;
            const label nbrNcmwPatchi =
                nbrCompleteMesh.boundary()[ncmwFvp.nbrPatchName()].index();

            checkNonConformalCoupledPatchOrdering
            (
                {-labelMax, labelMax},
                completeMesh.boundary()[ncmwPatchi],
                nbrCompleteMesh.boundary()[nbrNcmwPatchi],
                completeMesh.polyFacesBf()[ncmwPatchi],
                nbrCompleteMesh.polyFacesBf()[ncmwPatchi]
            );
        }

        // Error patches
        if (isA<nonConformalErrorFvPatch>(fvp))
        {
            const label ncePatchi = patchi;

            checkNonConformalErrorPatchOrdering
            (
                -labelMax,
                completeMesh.boundary()[ncePatchi],
                completeMesh.polyFacesBf()[ncePatchi]
            );
        }
    }
}


void checkProcMeshesOrdering
(
    const PtrList<fvMesh>& procMeshes,
    const multiDomainDecomposition& regionMeshes
)
{
    forAll(procMeshes, proci)
    {
        forAll(procMeshes[proci].boundary(), patchi)
        {
            const fvPatch& fvp =
                procMeshes[proci].boundary()[patchi];

            // Cyclic patches
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
                        procMeshes[nbrProci].boundary();

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
                    procMeshes[proci].boundary()[nccPatchi],
                    procMeshes[nbrProci].boundary()[nbrNccPatchi],
                    procMeshes[proci].polyFacesBf()[nccPatchi],
                    procMeshes[nbrProci].polyFacesBf()[nbrNccPatchi]
                );
            }

            // Mapped patches
            if
            (
                isA<nonConformalMappedWallFvPatch>(fvp)
             && refCast<const nonConformalMappedWallFvPatch>(fvp).owner()
            )
            {
                const label ncmwPatchi = patchi;

                const nonConformalMappedWallFvPatch& ncmwFvp =
                    refCast<const nonConformalMappedWallFvPatch>(fvp);

                typedef
                    nonConformalMappedPolyFacesFvsPatchLabelField
                    NcmpfFvsplf;

                const NcmpfFvsplf& polyFacesPf =
                    refCast<const NcmpfFvsplf>
                    (procMeshes[proci].polyFacesBf()[ncmwPatchi]);

                forAll(procMeshes, nbrProci)
                {
                    const fvMesh& nbrProcMesh =
                        regionMeshes[ncmwFvp.nbrRegionName()]()
                       .procMeshes()[nbrProci];

                    if (nbrProcMesh.conformal()) continue;

                    const label nbrNcmwPatchi =
                        nbrProcMesh
                       .boundary()[ncmwFvp.nbrPatchName()]
                       .index();

                    const nonConformalMappedWallFvPatch& nbrNcmwFvp =
                        refCast<const nonConformalMappedWallFvPatch>
                        (nbrProcMesh.boundary()[nbrNcmwPatchi]);

                    const NcmpfFvsplf& nbrPolyFacesPf =
                        refCast<const NcmpfFvsplf>
                        (nbrProcMesh.polyFacesBf()[nbrNcmwPatchi]);

                    checkNonConformalCoupledPatchOrdering
                    (
                        {proci, nbrProci},
                        ncmwFvp,
                        nbrNcmwFvp,
                        SubList<label>
                        (
                            procMeshes[proci].polyFacesBf()[ncmwPatchi],
                            polyFacesPf.procSizes()[nbrProci],
                            polyFacesPf.procOffsets()[nbrProci]
                        ),
                        SubList<label>
                        (
                            nbrProcMesh.polyFacesBf()[nbrNcmwPatchi],
                            nbrPolyFacesPf.procSizes()[proci],
                            nbrPolyFacesPf.procOffsets()[proci]
                        )
                    );
                }
            }

            // Error patches
            if (isA<nonConformalErrorFvPatch>(fvp))
            {
                const label ncePatchi = patchi;

                checkNonConformalErrorPatchOrdering
                (
                    proci,
                    procMeshes[proci].boundary()[ncePatchi],
                    procMeshes[proci].polyFacesBf()[ncePatchi]
                );
            }
        }
    }
}


}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::domainDecomposition::sortReconstructNonConformalCyclicAddressing_ =
    Foam::debug::optimisationSwitch
    (
        "sortReconstructNonConformalCyclicAddressing",
        0
    );


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


Foam::labelList Foam::domainDecomposition::completeFaceAddressing() const
{
    labelList result (completeMesh().nFaces(), -labelMax);

    forAll(procMeshes_, proci)
    {
        forAll(procFaceAddressing()[proci], procFacei)
        {
            const label facei =
                mag(procFaceAddressing()[proci][procFacei]) - 1;

            result[facei] = result[facei] == -labelMax ? procFacei : -1;
        }
    }

    return result;
}


Foam::List<Foam::domainDecomposition::labelPairLabelTable>
Foam::domainDecomposition::nonConformalCyclicProcCyclics() const
{
    List<labelPairLabelTable> result(completeMesh().boundary().size());

    forAll(procMeshes_, proci)
    {
        const fvMesh& procMesh = procMeshes_[proci];

        forAll(procMesh.boundary(), procPatchi)
        {
            const fvPatch& fvp = procMesh.boundary()[procPatchi];

            if (isA<nonConformalCyclicFvPatch>(fvp))
            {
                result[procPatchi].insert
                (
                    labelPair(proci, proci),
                    procPatchi
                );
            }

            if (isA<nonConformalProcessorCyclicFvPatch>(fvp))
            {
                const processorCyclicFvPatch& pcFvp =
                    refCast<const processorCyclicFvPatch>(fvp);

                result[pcFvp.referPatchIndex()].insert
                (
                    labelPair(proci, pcFvp.neighbProcNo()),
                    procPatchi
                );
            }
        }
    }

    return result;
}


Foam::PtrList<Foam::labelListList>
Foam::domainDecomposition::nonConformalMappedWallProcOffsets
(
    const bool appendSize
) const
{
    PtrList<labelListList> result(completeMesh().boundary().size());

    const surfaceLabelField::Boundary& polyFacesBf =
        completeMesh().polyFacesBf();

    forAll(completeMesh().boundary(), ncmwPatchi)
    {
        const fvPatch& fvp = completeMesh().boundary()[ncmwPatchi];

        if (!isA<nonConformalMappedWallFvPatch>(fvp)) continue;

        const nonConformalMappedWallFvPatch& ncmwFvp =
            refCast<const nonConformalMappedWallFvPatch>(fvp);

        const domainDecomposition& nbrDecomposition =
            regionMeshes_[ncmwFvp.nbrRegionName()]();

        const fvMesh& nbrCompleteMesh = nbrDecomposition.completeMesh();

        const label nbrNcmwPatchi =
            nbrCompleteMesh.boundary()[ncmwFvp.nbrPatchName()].index();

        const surfaceLabelField::Boundary& nbrPolyFacesBf =
            nbrCompleteMesh.polyFacesBf();

        result.set
        (
            ncmwPatchi,
            new labelListList(nProcs(), labelList(nProcs() + appendSize, 0))
        );

        // Determine the number of faces in each processor block of the
        // mapped patches
        labelListList& procNbrProcCounts = result[ncmwPatchi];
        forAll(polyFacesBf[ncmwPatchi], ncmwPatchFacei)
        {
            const label facei = polyFacesBf[ncmwPatchi][ncmwPatchFacei];
            const label celli = completeMesh().faceOwner()[facei];
            const label proci = cellProc_[celli];

            const label nbrFacei =
                nbrPolyFacesBf[nbrNcmwPatchi][ncmwPatchFacei];
            const label nbrCelli =
                nbrCompleteMesh.faceOwner()[nbrFacei];
            const label nbrProci = nbrDecomposition.cellProc()[nbrCelli];

            procNbrProcCounts[proci][nbrProci] ++;
        }

        // Convert the counts to cumulative sums (i.e., offsets)
        forAll(procNbrProcCounts, proci)
        {
            label count = 0;

            forAll(procNbrProcCounts[proci], nbrProci)
            {
                const label count0 = count;
                count += procNbrProcCounts[proci][nbrProci];
                procNbrProcCounts[proci][nbrProci] = count0;
            }

            if (appendSize)
            {
                procNbrProcCounts[proci].last() = count;
            }
        }
    }

    return result;
}


void Foam::domainDecomposition::decomposeNonConformalCyclicAddressing
(
    const label nccPatchi,
    const List<labelPairLabelTable>& nonConformalCyclicProcCyclics,
    List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
) const
{
    const surfaceLabelField::Boundary& polyFacesBf =
        completeMesh().polyFacesBf();

    const nonConformalCyclicFvPatch& nccFvp =
        refCast<const nonConformalCyclicFvPatch>
        (
            completeMesh().boundary()[nccPatchi]
        );

    if (!nccFvp.owner()) return;

    const label nbrNccPatchi = nccFvp.nbrPatchIndex();

    forAll(polyFacesBf[nccPatchi], nccPatchFacei)
    {
        const label facei = polyFacesBf[nccPatchi][nccPatchFacei];
        const label celli = completeMesh().faceOwner()[facei];
        const label proci = cellProc_[celli];

        const label nbrFacei = polyFacesBf[nbrNccPatchi][nccPatchFacei];
        const label nbrCelli = completeMesh().faceOwner()[nbrFacei];
        const label nbrProci = cellProc_[nbrCelli];

        const label procNccPatchi =
            nonConformalCyclicProcCyclics
            [nccPatchi][labelPair(proci, nbrProci)];
        const label nbrProcNccPatchi =
            nonConformalCyclicProcCyclics
            [nbrNccPatchi][labelPair(nbrProci, proci)];

        nonConformalProcFaceAddressingBf[proci][procNccPatchi]
            .append(nccPatchFacei + 1);
        nonConformalProcFaceAddressingBf[nbrProci][nbrProcNccPatchi]
            .append(nccPatchFacei + 1);
    }
}


void Foam::domainDecomposition::decomposeNonConformalMappedWallAddressing
(
    const label ncmwPatchi,
    PtrList<labelListList>& nonConformalMappedWallProcOffsets,
    List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
) const
{
    const surfaceLabelField::Boundary& polyFacesBf =
        completeMesh().polyFacesBf();

    const nonConformalMappedWallFvPatch& ncmwFvp =
        refCast<const nonConformalMappedWallFvPatch>
        (
            completeMesh().boundary()[ncmwPatchi]
        );

    const domainDecomposition& nbrDecomposition =
        regionMeshes_[ncmwFvp.nbrRegionName()]();

    const fvMesh& nbrCompleteMesh = nbrDecomposition.completeMesh();

    const surfaceLabelField::Boundary& nbrPolyFacesBf =
        nbrCompleteMesh.polyFacesBf();

    const label nbrNcmwPatchi =
        nbrCompleteMesh.boundary()[ncmwFvp.nbrPatchName()].index();

    // Resize the addressing as necessary
    labelListList& procOffsets =
        nonConformalMappedWallProcOffsets[ncmwPatchi];
    forAll(procOffsets, proci)
    {
        nonConformalProcFaceAddressingBf[proci][ncmwPatchi]
            .resize(procOffsets[proci].last());
    }

    // Insert the poly face addressing into the result. Use the
    // procOffsets array as the index into each processor block.
    forAll(polyFacesBf[ncmwPatchi], ncmwPatchFacei)
    {
        const label facei = polyFacesBf[ncmwPatchi][ncmwPatchFacei];
        const label celli = completeMesh().faceOwner()[facei];
        const label proci = cellProc_[celli];

        const label nbrFacei = nbrPolyFacesBf[nbrNcmwPatchi][ncmwPatchFacei];
        const label nbrCelli = nbrCompleteMesh.faceOwner()[nbrFacei];
        const label nbrProci = nbrDecomposition.cellProc()[nbrCelli];

        nonConformalProcFaceAddressingBf
            [proci][ncmwPatchi][procOffsets[proci][nbrProci] ++] =
            ncmwPatchFacei + 1;
    }
}


void Foam::domainDecomposition::decomposeNonConformalErrorAddressing
(
    const label ncePatchi,
    List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
) const
{
    const surfaceLabelField::Boundary& polyFacesBf =
        completeMesh().polyFacesBf();

    forAll(polyFacesBf[ncePatchi], ncePatchFacei)
    {
        const label facei = polyFacesBf[ncePatchi][ncePatchFacei];
        const label celli = completeMesh().faceOwner()[facei];
        const label proci = cellProc_[celli];

        nonConformalProcFaceAddressingBf[proci][ncePatchi]
            .append(ncePatchFacei + 1);
    }
}


void Foam::domainDecomposition::reconstructNonConformalCyclicAddressing
(
    const label nccPatchi,
    const List<labelPairLabelTable>& nonConformalCyclicProcCyclics,
    List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
) const
{
    const nonConformalCyclicFvPatch& nccFvp =
        refCast<const nonConformalCyclicFvPatch>
        (
            completeMesh().boundary()[nccPatchi]
        );

    if (!nccFvp.owner()) return;

    const label nbrNccPatchi = nccFvp.nbrPatchIndex();

    label nccPatchFacei = 0;
    labelPairLabelTable procNccPatchFaceis;
    forAllConstIter
    (
        labelPairLabelTable,
        nonConformalCyclicProcCyclics[nccPatchi],
        iter
    )
    {
        procNccPatchFaceis.insert(iter.key(), 0);
    }

    while (true)
    {
        labelPair procNbrProc(labelMax, labelMax);
        labelPair faceNbrFace(labelMax, labelMax);

        forAllConstIter(labelPairLabelTable, procNccPatchFaceis, iter)
        {
            const label proci = iter.key().first();
            const label nbrProci = iter.key().second();

            const labelPair procNbrProcStar(proci, nbrProci);
            const labelPair nbrProcProcStar(nbrProci, proci);

            const label procNccPatchi =
                nonConformalCyclicProcCyclics[nccPatchi][procNbrProcStar];
            const label nbrProcNccPatchi =
                nonConformalCyclicProcCyclics[nbrNccPatchi][nbrProcProcStar];

            const label procNccPatchFacei = iter();
            const label size =
                procMeshes_[proci].polyFacesBf()[procNccPatchi].size();

            if (procNccPatchFacei >= size) continue;

            const label procFacei =
                procMeshes_[proci].polyFacesBf()
                [procNccPatchi][procNccPatchFacei];
            const label nbrProcFacei =
                procMeshes_[nbrProci].polyFacesBf()
                [nbrProcNccPatchi][procNccPatchFacei];

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
                nonConformalCyclicProcCyclics[nccPatchi][procNbrProc];
            const label nbrProcNccPatchi =
                nonConformalCyclicProcCyclics[nbrNccPatchi][nbrProcProc];

            nonConformalProcFaceAddressingBf[proci][procNccPatchi]
                .append(nccPatchFacei + 1);
            nonConformalProcFaceAddressingBf[nbrProci][nbrProcNccPatchi]
                .append(nccPatchFacei + 1);

            nccPatchFacei ++;
            procNccPatchFaceis[procNbrProc] ++;
        }
    }
}


void Foam::domainDecomposition::sortReconstructNonConformalCyclicAddressing
(
    const label nccPatchi,
    const List<labelPairLabelTable>& nonConformalCyclicProcCyclics,
    List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
) const
{
    const nonConformalCyclicFvPatch& nccFvp =
        refCast<const nonConformalCyclicFvPatch>
        (
            completeMesh().boundary()[nccPatchi]
        );

    if (!nccFvp.owner()) return;

    const label nbrNccPatchi = nccFvp.nbrPatchIndex();

    // Resize the relevant patches
    forAllConstIter
    (
        labelPairLabelTable,
        nonConformalCyclicProcCyclics[nccPatchi],
        iter
    )
    {
        const label proci = iter.key().first();
        const label nbrProci = iter.key().second();

        const label procNccPatchi = iter();
        const label nbrProcNccPatchi =
            nonConformalCyclicProcCyclics
            [nbrNccPatchi][labelPair(nbrProci, proci)];

        const labelUList& procPolyFacesPf =
            procMeshes_[proci].polyFacesBf()[procNccPatchi];
        const labelUList& nbrProcPolyFacesPf =
            procMeshes_[nbrProci].polyFacesBf()[nbrProcNccPatchi];

        nonConformalProcFaceAddressingBf[proci][procNccPatchi]
            .resize(procPolyFacesPf.size());
        nonConformalProcFaceAddressingBf[nbrProci][nbrProcNccPatchi]
            .resize(nbrProcPolyFacesPf.size());
    }

    // Obtain references to the original patches
    const fvPatch& origFvp = nccFvp.origPatch();
    const fvPatch& nbrOrigFvp = nccFvp.nbrPatch().origPatch();

    // Create a list of "indices", containing every label relevant to each
    // non-conformal face
    DynamicList<FixedList<label, 7>> indices;
    forAllConstIter
    (
        labelPairLabelTable,
        nonConformalCyclicProcCyclics[nccPatchi],
        iter
    )
    {
        const label proci = iter.key().first();
        const label nbrProci = iter.key().second();

        const label procNccPatchi = iter();
        const label nbrProcNccPatchi =
            nonConformalCyclicProcCyclics
            [nbrNccPatchi][labelPair(nbrProci, proci)];

        const labelUList& procPolyFacesPf =
            procMeshes_[proci].polyFacesBf()[procNccPatchi];
        const labelUList& nbrProcPolyFacesPf =
            procMeshes_[nbrProci].polyFacesBf()[nbrProcNccPatchi];

        forAll(procPolyFacesPf, procPatchFacei)
        {
            const label procPolyFacei =
                procPolyFacesPf[procPatchFacei];
            const label nbrProcPolyFacei =
                nbrProcPolyFacesPf[procPatchFacei];

            const label completePolyFacei =
                procFaceAddressing_[proci][procPolyFacei] - 1;
            const label nbcCompletePolyFacei =
                procFaceAddressing_[nbrProci][nbrProcPolyFacei] - 1;

            indices.append
            ({
                proci,
                nbrProci,
                procNccPatchi,
                nbrProcNccPatchi,
                completePolyFacei - origFvp.start(),
                nbcCompletePolyFacei - nbrOrigFvp.start(),
                procPatchFacei
            });
        }
    }

    // Sort the indices by the owner poly face, then by the neighbour poly face
    Foam::stableSort
    (
        indices,
        [](const FixedList<label, 7>& a, const FixedList<label, 7>& b)
        {
            return labelPair(a[4], a[5]) < labelPair(b[4], b[5]);
        }
    );

    // Unpack into the addressing
    forAll(indices, i)
    {
        const label proci = indices[i][0];
        const label nbrProci = indices[i][1];
        const label procNccPatchi = indices[i][2];
        const label nbrProcNccPatchi = indices[i][3];
        //const label completePolyPatchFacei = indices[i][4];
        //const label nbrCompletePolyPatchFacei = indices[i][5];
        const label procPatchFacei = indices[i][6];

        nonConformalProcFaceAddressingBf
            [proci][procNccPatchi][procPatchFacei] = i + 1;
        nonConformalProcFaceAddressingBf
            [nbrProci][nbrProcNccPatchi][procPatchFacei] = i + 1;
    }
}


void Foam::domainDecomposition::reconstructNonConformalMappedWallAddressing
(
    const label ncmwPatchi,
    List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
) const
{
    const nonConformalMappedWallFvPatch& ncmwFvp =
        refCast<const nonConformalMappedWallFvPatch>
        (
            completeMesh().boundary()[ncmwPatchi]
        );

    const bool owner = ncmwFvp.owner();

    const domainDecomposition& nbrDecomposition =
        regionMeshes_[ncmwFvp.nbrRegionName()]();

    const label nbrNcmwPatchi =
        nbrDecomposition.completeMesh()
       .boundary()[ncmwFvp.nbrPatchName()]
       .index();

    auto calcProcOffsets = []
    (
        const domainDecomposition& meshes,
        const label ncmwPatchi
    )
    {
        labelListList result(meshes.nProcs());

        forAll(meshes.procMeshes(), proci)
        {
            typedef
                nonConformalMappedPolyFacesFvsPatchLabelField
                NcmpfFvsplf;

            const NcmpfFvsplf& polyFacesPf =
                refCast<const NcmpfFvsplf>
                (meshes.procMeshes()[proci].polyFacesBf()[ncmwPatchi]);

            result[proci] = polyFacesPf.procOffsets();
            result[proci].append(polyFacesPf.size());
        }

        return result;
    };
    const labelListList procOffsets =
        calcProcOffsets(*this, ncmwPatchi);
    const labelListList nbrProcOffsets =
        calcProcOffsets(nbrDecomposition, nbrNcmwPatchi);

    label ncmwPatchFacei = 0;
    labelListList procNcmwPatchFaceis(nProcs(), labelList(nProcs(), 0));

    while (true)
    {
        labelPair procNbrProc(labelMax, labelMax);
        labelPair faceNbrFace(labelMax, labelMax);

        forAll(procNcmwPatchFaceis, proci)
        {
            forAll(procNcmwPatchFaceis[proci], nbrProci)
            {
                const labelPair procNbrProcStar(proci, nbrProci);

                const label procNcmwPatchFacei =
                    procNcmwPatchFaceis[proci][nbrProci];
                const label size =
                    procOffsets[proci][nbrProci + 1]
                  - procOffsets[proci][nbrProci];

                if (procNcmwPatchFacei >= size) continue;

                const label procFacei =
                    procMeshes_[proci].polyFacesBf()
                    [ncmwPatchi]
                    [procNcmwPatchFacei + procOffsets[proci][nbrProci]];
                const label nbrProcFacei =
                    nbrDecomposition.procMeshes()[nbrProci].polyFacesBf()
                    [nbrNcmwPatchi]
                    [procNcmwPatchFacei + nbrProcOffsets[nbrProci][proci]];

                const labelPair faceNbrFaceStar
                (
                    procFaceAddressing_[proci][procFacei] - 1,
                    nbrDecomposition
                   .procFaceAddressing_[nbrProci][nbrProcFacei] - 1
                );

                if
                (
                    owner
                  ? faceNbrFace > faceNbrFaceStar
                  : reverse(faceNbrFace) > reverse(faceNbrFaceStar)
                )
                {
                    procNbrProc = procNbrProcStar;
                    faceNbrFace = faceNbrFaceStar;
                }
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

            nonConformalProcFaceAddressingBf[proci][ncmwPatchi]
                .append(ncmwPatchFacei + 1);

            ncmwPatchFacei ++;
            procNcmwPatchFaceis[proci][nbrProci] ++;
        }
    }
}


void Foam::domainDecomposition::reconstructNonConformalErrorAddressing
(
    const label ncePatchi,
    List<List<DynamicList<label>>>& nonConformalProcFaceAddressingBf
) const
{
    label ncePatchFacei = 0;
    labelList procNcePatchFaceis(nProcs(), 0);

    while (true)
    {
        label facei = labelMax, proci = labelMax;

        forAll(procNcePatchFaceis, procStari)
        {
            const label size =
                procMeshes_[procStari].polyFacesBf()[ncePatchi].size();

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
            nonConformalProcFaceAddressingBf[proci][ncePatchi]
                .append(ncePatchFacei + 1);

            ncePatchFacei ++;
            procNcePatchFaceis[proci] ++;
        }
    }
}


Foam::List<Foam::List<Foam::DynamicList<Foam::label>>>
Foam::domainDecomposition::nonConformalProcFaceAddressingBf() const
{
    validateComplete();
    validateProcs();

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
        const List<labelPairLabelTable> nonConformalCyclicProcCyclics =
            this->nonConformalCyclicProcCyclics();

        PtrList<labelListList> nonConformalMappedWallProcOffsets =
            this->nonConformalMappedWallProcOffsets(true);

        // Decompose non-conformal addressing
        forAll(completeMesh().boundary(), ncPatchi)
        {
            const fvPatch& fvp = completeMesh().boundary()[ncPatchi];

            if (isA<nonConformalCyclicFvPatch>(fvp))
            {
                decomposeNonConformalCyclicAddressing
                (
                    ncPatchi,
                    nonConformalCyclicProcCyclics,
                    result
                );
            }

            if (isA<nonConformalMappedWallFvPatch>(fvp))
            {
                decomposeNonConformalMappedWallAddressing
                (
                    ncPatchi,
                    nonConformalMappedWallProcOffsets,
                    result
                );
            }

            if (isA<nonConformalErrorFvPatch>(fvp))
            {
                decomposeNonConformalErrorAddressing(ncPatchi, result);
            }
        }
    }
    else // if (!procsConformal())
    {
        const List<labelPairLabelTable> nonConformalCyclicProcCyclics =
            this->nonConformalCyclicProcCyclics();

        // Reconstruct non-conformal addressing
        forAll(completeMesh().boundary(), ncPatchi)
        {
            const fvPatch& fvp = completeMesh().boundary()[ncPatchi];

            if (isA<nonConformalCyclicFvPatch>(fvp))
            {
                if (sortReconstructNonConformalCyclicAddressing_)
                {
                    sortReconstructNonConformalCyclicAddressing
                    (
                        ncPatchi,
                        nonConformalCyclicProcCyclics,
                        result
                    );
                }
                else
                {
                    reconstructNonConformalCyclicAddressing
                    (
                        ncPatchi,
                        nonConformalCyclicProcCyclics,
                        result
                    );
                }
            }

            if (isA<nonConformalMappedWallFvPatch>(fvp))
            {
                reconstructNonConformalMappedWallAddressing(ncPatchi, result);
            }

            if (isA<nonConformalErrorFvPatch>(fvp))
            {
                reconstructNonConformalErrorAddressing(ncPatchi, result);
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

        forAll(procMesh.boundary(), procNcPatchi)
        {
            const fvPatch& fvp = procMesh.boundary()[procNcPatchi];

            if (!isA<nonConformalFvPatch>(fvp)) continue;

            const label completeNcPatchi =
                isA<processorCyclicFvPatch>(fvp)
              ? refCast<const processorCyclicFvPatch>(fvp)
               .referPatchIndex()
              : procNcPatchi;

            const label size =
                max
                (
                    max(mag(faceAddressingBf[procNcPatchi])),
                    polyFacesBf[completeNcPatchi].size()
                );

            polyFacesBf[completeNcPatchi].resize(size, -1);
            polyFacesBf[completeNcPatchi].labelField::rmap
            (
                mag
                (
                    labelField
                    (
                        procFaceAddressing_[proci],
                        procMesh.polyFacesBf()[procNcPatchi]
                    )
                ) - 1,
                mag(faceAddressingBf[procNcPatchi]) - 1
            );

            // Set dummy data for the face geometry. This should not be
            // used during reconstruction.
            Sf.boundaryFieldRef()[completeNcPatchi].resize(size, Zero);
            Cf.boundaryFieldRef()[completeNcPatchi].resize(size, Zero);
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

    completeMesh_->setPolyFacesBfInstance(procMeshes_[0].polyFacesBfInstance());

    if (debug)
    {
        checkCompleteMeshOrdering(completeMesh(), regionMeshes_);
    }
}


void Foam::domainDecomposition::unconformProcs()
{
    const labelList completeFaceAddressing =
        this->completeFaceAddressing();

    const PtrList<labelListList> nonConformalMappedWallProcOffsets =
        this->nonConformalMappedWallProcOffsets(false);

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

        forAll(procMesh.boundary(), procNcPatchi)
        {
            const fvPatch& fvp = procMesh.boundary()[procNcPatchi];

            if (!isA<nonConformalFvPatch>(fvp)) continue;

            const label completeNcPatchi =
                isA<processorCyclicFvPatch>(fvp)
              ? refCast<const processorCyclicFvPatch>(fvp)
               .referPatchIndex()
              : procNcPatchi;

            polyFacesBf[procNcPatchi] =
                labelField
                (
                    completeFaceAddressing,
                    labelField
                    (
                        completeMesh().polyFacesBf()[completeNcPatchi],
                        mag(faceAddressingBf[procNcPatchi]) - 1
                    )
                );

            if (isA<nonConformalMappedWallFvPatch>(fvp))
            {
                typedef
                    nonConformalMappedPolyFacesFvsPatchLabelField
                    NcmpfFvsplf;

                refCast<NcmpfFvsplf>(polyFacesBf[procNcPatchi]).procOffsets() =
                    nonConformalMappedWallProcOffsets[procNcPatchi][proci];
            }

            const label size = polyFacesBf[procNcPatchi].size();

            // Set dummy data for the face geometry. This should not be
            // used during decomposition.
            Sf.boundaryFieldRef()[procNcPatchi].resize(size, Zero);
            Cf.boundaryFieldRef()[procNcPatchi].resize(size, Zero);
        }

        procMesh.unconform
        (
            polyFacesBf,
            Sf,
            Cf,
            NullObjectRef<surfaceScalarField>(),
            false
        );

        procMesh.setPolyFacesBfInstance(completeMesh().polyFacesBfInstance());
    }

    // Check ordering
    if (debug)
    {
        checkProcMeshesOrdering(procMeshes(), regionMeshes_);
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
