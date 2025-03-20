/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "fvMeshStitcher.H"
#include "globalIndex.H"
#include "fvcSurfaceIntegrate.H"
#include "MultiRegionList.H"
#include "meshObjects.H"
#include "movingWallVelocityFvPatchVectorField.H"
#include "movingWallSlipVelocityFvPatchVectorField.H"
#include "nonConformalBoundary.H"
#include "nonConformalCyclicFvPatch.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "nonConformalErrorFvPatch.H"
#include "nonConformalCyclicFvPatch.H"
#include "nonConformalMappedWallFvPatch.H"
#include "nonConformalMappedPolyFacesFvsPatchLabelField.H"
#include "nonConformalProcessorCyclicFvPatch.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "syncTools.H"
#include "surfaceInterpolate.H"
#include "surfaceToVolVelocity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class FaceList, class PointField>
edge meshEdge
(
    const PrimitivePatch<FaceList, PointField>& p,
    const label edgei
)
{
    return edge
    (
        p.meshPoints()[p.edges()[edgei][0]],
        p.meshPoints()[p.edges()[edgei][1]]
    );
}


bool any(const boolList& l)
{
    bool result = false;
    forAll(l, i)
    {
        if (l[i])
        {
            result = true;
            break;
        }
    }
    return returnReduce(result, andOp<bool>());
}

}


// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

const Foam::scalar Foam::fvMeshStitcher::minWarnProjectedVolumeFraction_ = 0.3;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshStitcher, 0);
    defineRunTimeSelectionTable(fvMeshStitcher, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshStitcher::regionNames(const fvMesh& mesh, wordHashSet& names)
{
    names.insert(mesh.name());

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];

        if (!isA<nonConformalMappedWallFvPatch>(fvp)) continue;

        const nonConformalMappedWallFvPatch& ncmwFvp =
            refCast<const nonConformalMappedWallFvPatch>(fvp);

        if (!names.found(ncmwFvp.nbrRegionName()))
        {
            regionNames(ncmwFvp.nbrMesh(), names);
        }
    }
}


Foam::wordList Foam::fvMeshStitcher::regionNames() const
{
    wordHashSet names;
    regionNames(mesh_, names);
    return names.toc();
}


Foam::UPtrList<Foam::fvMesh> Foam::fvMeshStitcher::regionMeshes()
{
    const wordList names(regionNames());

    UPtrList<fvMesh> result(names.size());
    forAll(names, i)
    {
        result.set
        (
            i,
            &mesh_.time().lookupObjectRef<fvMesh>(names[i])
        );
    }

    return result;
}


Foam::UPtrList<const Foam::fvMesh> Foam::fvMeshStitcher::regionMeshes() const
{
    const wordList names(regionNames());

    UPtrList<const fvMesh> result(names.size());
    forAll(names, i)
    {
        result.set
        (
            i,
            &mesh_.time().lookupObject<fvMesh>(names[i])
        );
    }

    return result;
}


Foam::IOobject Foam::fvMeshStitcher::polyFacesBfIO(const fvMesh& mesh)
{
    return IOobject
    (
        "polyFaces",
        mesh.time().findInstance
        (
            mesh.dbDir()/fvMesh::typeName,
            "polyFaces",
            IOobject::READ_IF_PRESENT
        ),
        fvMesh::typeName,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );
}


bool Foam::fvMeshStitcher::loadPolyFacesBf
(
    IOobject& polyFacesBfIO,
    SurfaceFieldBoundary<label>& polyFacesBf
)
{
    bool loaded = false;

    // If the poly faces are in the region cache then use them and remove them.
    // This mesh is now available, so we can look it up to access its poly
    // faces from now on.
    auto polyFacesBfIOIter = regionPolyFacesBfIOs_.find(mesh_.name());
    auto polyFacesBfIter = regionPolyFacesBfs_.find(mesh_.name());
    if (polyFacesBfIOIter != regionPolyFacesBfIOs_.end())
    {
        polyFacesBfIO =
            autoPtr<IOobject>
            (
                regionPolyFacesBfIOs_.remove(polyFacesBfIOIter)
            )();

        polyFacesBf.reset
        (
            tmp<SurfaceFieldBoundary<label>>
            (
                regionPolyFacesBfs_.remove(polyFacesBfIter)
            )
        );

        loaded = true;
    }

    // If the faces are not in the region cache, then try and read them
    if (!loaded)
    {
        typeIOobject<surfaceLabelField> io
        (
            fvMeshStitcher::polyFacesBfIO(mesh_)
        );

        if (io.headerOk())
        {
            polyFacesBfIO = io;

            polyFacesBf.reset
            (
                surfaceLabelField(polyFacesBfIO, mesh_).boundaryField()
            );

            loaded = true;
        }
    }

    // If nothing is available then return
    if (!loaded) return false;

    // If the poly faces are available, then make sure all connected regions'
    // poly faces have been read and are in the cache
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!isA<nonConformalMappedWallFvPatch>(fvp)) continue;

        const nonConformalMappedWallFvPatch& ncmwFvp =
            refCast<const nonConformalMappedWallFvPatch>(fvp);

        if (!ncmwFvp.haveNbr()) continue;

        if (regionPolyFacesBfs_.found(ncmwFvp.nbrMesh().name())) continue;

        IOobject io = fvMeshStitcher::polyFacesBfIO(ncmwFvp.nbrMesh());

        regionPolyFacesBfIOs_.insert
        (
            ncmwFvp.nbrMesh().name(),
            new IOobject(io)
        );

        regionPolyFacesBfs_.insert
        (
            ncmwFvp.nbrMesh().name(),
            new surfaceLabelField::Boundary
            (
                surfaceLabelField::null(),
                surfaceLabelField(io, ncmwFvp.nbrMesh()).boundaryField()
            )
        );
    }

    return true;
}


const Foam::surfaceLabelField::Boundary& Foam::fvMeshStitcher::getPolyFacesBf
(
    const word& regionName
) const
{
    if (regionPolyFacesBfs_.found(regionName))
    {
        return *regionPolyFacesBfs_[regionName];
    }
    else
    {
        return mesh_.time().lookupObject<fvMesh>(regionName).polyFacesBf();
    }
}


void Foam::fvMeshStitcher::getOrigNbrBfs
(
    const SurfaceFieldBoundary<label>& polyFacesBf,
    const SurfaceFieldBoundary<vector>& SfBf,
    const SurfaceFieldBoundary<vector>& CfBf,
    tmp<SurfaceFieldBoundary<label>>& tOrigFacesNbrBf,
    tmp<SurfaceFieldBoundary<vector>>& tOrigSfNbrBf,
    tmp<SurfaceFieldBoundary<point>>& tOrigCfNbrBf
) const
{
    // Temporary full-sized surface fields for orig properties
    surfaceLabelField origFaces
    (
        surfaceLabelField::New
        (
            "origFaces",
            mesh_,
            dimensioned<label>(dimless, -1)
        )
    );
    surfaceVectorField origSf
    (
        surfaceVectorField::New
        (
            "origSf",
            mesh_,
            dimensionedVector("NaN", dimArea, vector::uniform(NaN))
        )
    );
    surfaceVectorField origCf
    (
        surfaceVectorField::New
        (
            "origCf",
            mesh_,
            dimensionedVector("NaN", dimLength, vector::uniform(NaN))
        )
    );

    // Communicate orig properties across cyclics and processor cyclics
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (isA<nonConformalCoupledFvPatch>(fvp))
        {
            const nonConformalCoupledFvPatch& nccFvp =
                refCast<const nonConformalCoupledFvPatch>(fvp);

            origFaces.boundaryFieldRef()[patchi] =
                polyFacesBf[patchi] - nccFvp.origPatch().start();
        }

        // (set Sf and Cf on all coupled patches, so that transform operations
        // in boundaryNeighbourField don't trigger a sigFpe)
        if (isA<coupledFvPatch>(fvp))
        {
            origSf.boundaryFieldRef()[patchi] =
                vectorField(mesh_.faceAreas(), polyFacesBf[patchi]);
            origCf.boundaryFieldRef()[patchi] =
                pointField(mesh_.faceCentres(), polyFacesBf[patchi]);
        }
    }
    origFaces.boundaryFieldRef() =
        origFaces.boundaryField().boundaryNeighbourField();
    origSf.boundaryFieldRef() =
        origSf.boundaryField().boundaryNeighbourField();
    origCf.boundaryFieldRef() =
        origCf.boundaryField().boundaryNeighbourField();

    // Communicate orig properties across mapped walls
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!isA<nonConformalMappedWallFvPatch>(fvp)) continue;

        const nonConformalMappedWallFvPatch& ncmwFvp =
            refCast<const nonConformalMappedWallFvPatch>(fvp);

        if (!ncmwFvp.haveNbr()) continue;

        const fvMesh& nbrMesh = ncmwFvp.nbrMesh();
        const nonConformalMappedWallFvPatch& nbrPatch = ncmwFvp.nbrPatch();

        const fvsPatchLabelField& nbrPolyFacesPf =
            getPolyFacesBf(nbrMesh.name())[nbrPatch.index()];

        origFaces.boundaryFieldRef()[patchi] =
            nonConformalMappedFvPatchBase::map
            (
                nbrPolyFacesPf,
                nbrPolyFacesPf - nbrPatch.origPatch().start(),
                polyFacesBf[patchi]
            );
        origSf.boundaryFieldRef()[patchi] =
            nonConformalMappedFvPatchBase::map
            (
                nbrPolyFacesPf,
                vectorField(nbrMesh.faceAreas(), nbrPolyFacesPf),
                polyFacesBf[patchi]
            );
        origCf.boundaryFieldRef()[patchi] =
            nonConformalMappedFvPatchBase::map
            (
                nbrPolyFacesPf,
                vectorField(nbrMesh.faceCentres(), nbrPolyFacesPf),
                polyFacesBf[patchi]
            );
    }

    // Correct transformed face centres. See note in fvMesh::unconform
    // regarding the transformation of cell centres. The same applies here.
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        fvsPatchVectorField& Cfp = origCf.boundaryFieldRef()[patchi];

        if (isA<nonConformalCoupledFvPatch>(fvp))
        {
            const nonConformalCoupledFvPatch& nccFvp =
                refCast<const nonConformalCoupledFvPatch>(fvp);

            nccFvp.transform().invTransform(Cfp, Cfp);
            nccFvp.transform().transformPosition(Cfp, Cfp);
        }

        if (isA<nonConformalMappedWallFvPatch>(fvp))
        {
            const nonConformalMappedWallFvPatch& ncmwFvp =
                refCast<const nonConformalMappedWallFvPatch>(fvp);

            if (ncmwFvp.haveNbr())
            {
                ncmwFvp.transform().invTransform(Cfp, Cfp);
                ncmwFvp.transform().transformPosition(Cfp, Cfp);
            }
        }
    }

    // Set the tmp boundary fields and discard the unused internal field
    tOrigFacesNbrBf =
        new surfaceLabelField::Boundary
        (
            surfaceLabelField::null(),
            origFaces.boundaryField()
        );
    tOrigSfNbrBf =
        new surfaceVectorField::Boundary
        (
            surfaceVectorField::null(),
            origSf.boundaryField()
        );
    tOrigCfNbrBf =
        new surfaceVectorField::Boundary
        (
            surfaceVectorField::null(),
            origCf.boundaryField()
        );
}


Foam::List<Foam::List<Foam::FixedList<Foam::label, 3>>>
Foam::fvMeshStitcher::procFacesToIndices
(
    const List<List<remote>>& faceOtherProcFaces,
    const bool owner
)
{
    // Compute interface sizes
    labelList is(Pstream::nProcs(), 0);
    forAll(faceOtherProcFaces, facei)
    {
        forAll(faceOtherProcFaces[facei], i)
        {
            const remote& otherProcFacei = faceOtherProcFaces[facei][i];

            is[otherProcFacei.proci] ++;
        }
    }

    // Allocate the indices
    List<List<FixedList<label, 3>>> indices(Pstream::nProcs());
    forAll(indices, proci)
    {
        indices[proci].resize(is[proci]);
    }

    // Set the indices
    is = 0;
    forAll(faceOtherProcFaces, facei)
    {
        forAll(faceOtherProcFaces[facei], i)
        {
            const remote& otherProcFacei = faceOtherProcFaces[facei][i];

            FixedList<label, 3>& index =
                indices[otherProcFacei.proci][is[otherProcFacei.proci] ++];

            index = {facei, otherProcFacei.elementi, i + 1};

            if (!owner) Swap(index[0], index[1]);
        }
    }

    // Sort to ensure ordering
    forAll(indices, proci)
    {
        sort(indices[proci]);
    }

    return indices;
}


void Foam::fvMeshStitcher::matchIndices
(
    const surfaceLabelField::Boundary& polyFacesBf,
    const surfaceLabelField::Boundary& origFacesNbrBf,
    List<List<FixedList<label, 3>>>& indices,
    const fvPatch& ncFvp,
    const polyPatch& origPp,
    const labelList& patchis,
    const labelList& patchOffsets,
    const labelList& patchSizes,
    const bool owner
)
{
    // Create what the indices should be
    List<List<FixedList<label, 3>>> indicesRef(indices.size());
    forAll(indices, proci)
    {
        const label patchi = patchis[proci];

        if (patchi != -1)
        {
            indicesRef[proci].resize(patchSizes[proci]);

            for (label indexi = 0; indexi < patchSizes[proci]; ++ indexi)
            {
                FixedList<label, 3>& indexRef = indicesRef[proci][indexi];

                const label patchFacei = indexi + patchOffsets[proci];

                indexRef =
                    {
                        polyFacesBf[patchi][patchFacei] - origPp.start(),
                        origFacesNbrBf[patchi][patchFacei],
                        0
                    };

                if (!owner) Swap(indexRef[0], indexRef[1]);
            }
        }
    }

    // Populate the indices with the index of the coupling
    label nCouplesRemoved = 0, nCouplesAdded = 0;
    forAll(indices, proci)
    {
        const label patchi = patchis[proci];

        if (patchi != -1)
        {
            label refi = 0, i = 0;

            DynamicList<FixedList<label, 3>> removedIndices;

            while
            (
                refi < indicesRef[proci].size()
             && i < indices[proci].size()
            )
            {
                const FixedList<label, 3> index
                ({
                    indices[proci][i][0],
                    indices[proci][i][1],
                    0
                });

                FixedList<label, 3>& indexRef = indicesRef[proci][refi];

                if (index < indexRef)
                {
                    nCouplesRemoved ++;
                    removedIndices.append
                    ({
                        indices[proci][i][0],
                        indices[proci][i][1],
                      - indices[proci][i][2]
                    });
                    i ++;
                }
                else if (index == indexRef)
                {
                    indexRef[2] = indices[proci][i][2];
                    refi ++;
                    i ++;
                }
                else // (index > indexRef)
                {
                    nCouplesAdded ++;
                    refi ++;
                }
            }

            nCouplesRemoved += min(indices[proci].size() - i, 0);
            nCouplesAdded += min(indicesRef[proci].size() - refi, 0);

            indicesRef[proci].append(removedIndices);
        }
    }

    // Report if changes have been made
    reduce(nCouplesRemoved, sumOp<label>());
    reduce(nCouplesAdded, sumOp<label>());
    if ((nCouplesRemoved || nCouplesAdded) && owner)
    {
        Info<< indent << nCouplesRemoved << '/' << nCouplesAdded
            << " small couplings removed/added to " << ncFvp.name()
            << endl;
    }

    // Set the indices to the correct values
    Swap(indices, indicesRef);
}


Foam::label Foam::fvMeshStitcher::nValidIndices
(
    const List<FixedList<label, 3>>& indices
)
{
    label n = indices.size();
    while (n > 0 && indices[n - 1][2] < 0) n --;
    return n;
}


void Foam::fvMeshStitcher::createCouplings
(
    surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField::Boundary& SfBf,
    surfaceVectorField::Boundary& CfBf,
    const tmp<surfaceVectorField::Boundary>& tOrigSfNbrBf,
    const tmp<surfaceVectorField::Boundary>& tOrigCfNbrBf,
    const List<List<FixedList<label, 3>>>& indices,
    const List<DynamicList<couple>>& couples,
    const polyPatch& origPp,
    const labelList& patchis,
    const labelList& patchOffsets,
    const bool owner
)
{
    // Add coupled faces into the non-conformal patches
    forAll(patchis, proci)
    {
        const label patchi = patchis[proci];
        const label patchOffset = patchOffsets[proci];

        if (patchi != -1)
        {
            forAll(indices[proci], indexi)
            {
                const label origFacei = indices[proci][indexi][!owner];
                const label i = indices[proci][indexi][2];

                const label patchFacei = i >= 0 ? indexi + patchOffset : -1;

                couple c;
                if (i != 0)
                {
                    c = couples[origFacei][mag(i) - 1];
                }
                else
                {
                    c =
                        couple
                        (
                            part
                            (
                                small*origPp.faceAreas()[origFacei],
                                origPp.faceCentres()[origFacei]
                            ),
                            part
                            (
                                small*tOrigSfNbrBf()[patchi][patchFacei],
                                tOrigCfNbrBf()[patchi][patchFacei]
                            )
                        );
                }

                // The two parts of the coupling. The projection is to the
                // neighbour, so the other-side is always taken from the
                // neighbouring patch faces.
                const part& pThis = c, pOther = owner ? c.nbr : c;

                // Remove the area from the corresponding original face
                if (i >= 0 || owner)
                {
                    part origP
                    (
                        SfBf[origPp.index()][origFacei],
                        CfBf[origPp.index()][origFacei]
                    );
                    origP -= pThis;
                    if (i < 0 && owner) origP += pOther;

                    SfBf[origPp.index()][origFacei] = origP.area;
                    CfBf[origPp.index()][origFacei] = origP.centre;
                }

                // Add the new coupled face
                if (i >= 0)
                {
                    polyFacesBf[patchi][patchFacei] =
                        origFacei + origPp.start();
                    SfBf[patchi][patchFacei] = pOther.area;
                    CfBf[patchi][patchFacei] = pOther.centre;
                }
            }
        }
    }
}


void Foam::fvMeshStitcher::createErrorAndEdgeParts
(
    surfaceVectorField::Boundary& SfBf,
    surfaceVectorField::Boundary& CfBf,
    List<part>& origEdgeParts,
    const patchToPatches::intersection& intersection,
    const polyPatch& origPp
)
{
    // Add error geometry to the original patches
    forAll(intersection.srcErrorParts(), origFacei)
    {
        part origP
        (
            SfBf[origPp.index()][origFacei],
            CfBf[origPp.index()][origFacei]
        );
        origP += intersection.srcErrorParts()[origFacei];

        SfBf[origPp.index()][origFacei] = origP.area;
        CfBf[origPp.index()][origFacei] = origP.centre;
    }

    // Store the orig edge parts
    forAll(intersection.srcEdgeParts(), origEdgei)
    {
        origEdgeParts[origEdgei] += intersection.srcEdgeParts()[origEdgei];
    }
}


void Foam::fvMeshStitcher::intersectNonConformalCyclic
(
    const nonConformalCyclicFvPatch& nccFvp,
    surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField::Boundary& SfBf,
    surfaceVectorField::Boundary& CfBf,
    const tmp<surfaceLabelField::Boundary>& tOrigFacesNbrBf,
    const tmp<surfaceVectorField::Boundary>& tOrigSfNbrBf,
    const tmp<surfaceVectorField::Boundary>& tOrigCfNbrBf,
    List<part>& origEdgeParts
) const
{
    // Alias the original poly patches
    const polyPatch& origPp = nccFvp.origPatch().patch();

    // Get the indices of the related (i.e., cyclic and processorCyclic)
    // non-conformal patches. Index based on the connected processor.
    labelList patchis(Pstream::nProcs(), -1);
    patchis[Pstream::myProcNo()] = nccFvp.index();
    forAll(mesh_.boundary(), patchj)
    {
        const fvPatch& fvp = mesh_.boundary()[patchj];

        if (isA<nonConformalProcessorCyclicFvPatch>(fvp))
        {
            const nonConformalProcessorCyclicFvPatch& ncpcFvp =
                refCast<const nonConformalProcessorCyclicFvPatch>(fvp);

            if (ncpcFvp.referPatchIndex() == nccFvp.index())
            {
                patchis[ncpcFvp.neighbProcNo()] = patchj;
            }
        }
    }

    // Get the intersection geometry
    const patchToPatches::intersection& intersection =
        nccFvp.owner()
      ? nccFvp.nonConformalCyclicPatch().intersection()
      : nccFvp.nbrPatch().nonConformalCyclicPatch().intersection();

    // Unpack the patchToPatch addressing into a list of indices
    List<List<FixedList<label, 3>>> indices =
        procFacesToIndices
        (
            nccFvp.owner()
          ? intersection.srcTgtProcFaces()
          : intersection.tgtSrcProcFaces(),
            nccFvp.owner()
        );

    // If addressing has been provided, then modify the indices to match
    if (tOrigFacesNbrBf.valid())
    {
        labelList patchSizes(Pstream::nProcs(), -1);
        forAll(patchis, proci)
        {
            if (patchis[proci] != -1)
            {
                patchSizes[proci] = polyFacesBf[patchis[proci]].size();
            }
        }

        matchIndices
        (
            polyFacesBf,
            tOrigFacesNbrBf(),
            indices,
            nccFvp,
            origPp,
            patchis,
            labelList(Pstream::nProcs(), 0),
            patchSizes,
            nccFvp.owner()
        );
    }

    // Check the existence of the non-conformal patches to be added into and
    // resize them as necessary
    forAll(patchis, proci)
    {
        const label patchi = patchis[proci];
        const label patchSize = nValidIndices(indices[proci]);

        if (patchi == -1 && patchSize)
        {
            FatalErrorInFunction
                << "Intersection generated coupling through "
                << nccFvp.name() << " between processes "
                << Pstream::myProcNo() << " and " << proci
                << " but a corresponding processor interface was not found"
                << exit(FatalError);
        }

        if (patchi != -1)
        {
            polyFacesBf[patchi].resize(patchSize);
            SfBf[patchi].resize(patchSize);
            CfBf[patchi].resize(patchSize);
        }
    }

    // Create couplings by transferring geometry from the original to the
    // non-conformal patches
    createCouplings
    (
        polyFacesBf,
        SfBf,
        CfBf,
        tOrigSfNbrBf,
        tOrigCfNbrBf,
        indices,
        nccFvp.owner() ? intersection.srcCouples() : intersection.tgtCouples(),
        origPp,
        patchis,
        labelList(Pstream::nProcs(), 0),
        nccFvp.owner()
    );

    // Add error geometry to the original patches and store the edge parts
    if (nccFvp.owner())
    {
        createErrorAndEdgeParts
        (
            SfBf,
            CfBf,
            origEdgeParts,
            intersection,
            origPp
        );
    }
}


void Foam::fvMeshStitcher::intersectNonConformalMappedWall
(
    const nonConformalMappedWallFvPatch& ncmwFvp,
    surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField::Boundary& SfBf,
    surfaceVectorField::Boundary& CfBf,
    const tmp<surfaceLabelField::Boundary>& tOrigFacesNbrBf,
    const tmp<surfaceVectorField::Boundary>& tOrigSfNbrBf,
    const tmp<surfaceVectorField::Boundary>& tOrigCfNbrBf,
    List<part>& origEdgeParts
) const
{
    // Alias the original poly patch
    const polyPatch& origPp = ncmwFvp.origPatch().patch();

    // Get the intersection geometry
    const patchToPatches::intersection& intersection =
        ncmwFvp.owner()
      ? ncmwFvp.nonConformalMappedWallPatch().intersection()
      : ncmwFvp.nbrPatch().nonConformalMappedWallPatch().intersection();

    // Unpack the patchToPatch addressing into a list of indices
    List<List<FixedList<label, 3>>> indices =
        procFacesToIndices
        (
            ncmwFvp.owner()
          ? intersection.srcTgtProcFaces()
          : intersection.tgtSrcProcFaces(),
            ncmwFvp.owner()
        );

    // Obtain the label patch field for the mapped wall poly faces
    nonConformalMappedPolyFacesFvsPatchLabelField& polyFacesPf =
        refCast<nonConformalMappedPolyFacesFvsPatchLabelField>
        (
            polyFacesBf[ncmwFvp.index()]
        );

    // If addressing has been provided, then modify the indices to match
    if (tOrigFacesNbrBf.valid())
    {
        matchIndices
        (
            polyFacesBf,
            tOrigFacesNbrBf(),
            indices,
            ncmwFvp,
            origPp,
            labelList(Pstream::nProcs(), ncmwFvp.index()),
            polyFacesPf.procOffsets(),
            polyFacesPf.procSizes(),
            ncmwFvp.owner()
        );
    }

    // Set the processor offsets within the mapped polyFacesBf patch field
    {
        label count = 0;

        forAll(indices, proci)
        {
            polyFacesPf.procOffsets()[proci] = count;

            count += nValidIndices(indices[proci]);
        }

        polyFacesPf.resize(count);
        SfBf[polyFacesPf.patch().index()].resize(count);
        CfBf[polyFacesPf.patch().index()].resize(count);
    };

    // Create a coupling by transferring geometry from the original to the
    // non-conformal patch
    createCouplings
    (
        polyFacesBf,
        SfBf,
        CfBf,
        tOrigSfNbrBf,
        tOrigCfNbrBf,
        indices,
        ncmwFvp.owner() ? intersection.srcCouples() : intersection.tgtCouples(),
        origPp,
        labelList(Pstream::nProcs(), ncmwFvp.index()),
        polyFacesPf.procOffsets(),
        ncmwFvp.owner()
    );

    // Add error geometry to the original patch and store the edge parts
    if (ncmwFvp.owner())
    {
        createErrorAndEdgeParts
        (
            SfBf,
            CfBf,
            origEdgeParts,
            intersection,
            origPp
        );
    }
}


Foam::List<Foam::fvMeshStitcher::part>
Foam::fvMeshStitcher::calculateOwnerOrigBoundaryEdgeParts
(
    const List<List<part>>& patchEdgeParts
) const
{
    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList& ownerOrigBoundaryPointMeshPoint =
        ncb.ownerOrigBoundaryPointMeshPoint();
    const labelList& ownerOrigBoundaryEdgeMeshEdge =
        ncb.ownerOrigBoundaryEdgeMeshEdge();
    const edgeList& ownerOrigBoundaryEdges = ncb.ownerOrigBoundaryEdges();
    const edgeList& ownerOrigBoundaryMeshEdges =
        ncb.ownerOrigBoundaryMeshEdges();

    // Sum up boundary edge parts, being careful with signs
    labelList ownerOrigBoundaryEdgeNParts(ownerOrigBoundaryEdges.size(), 0);
    List<part> ownerOrigBoundaryEdgeParts(ownerOrigBoundaryEdges.size());
    forAll(patchEdgeParts, patchi)
    {
        if (patchEdgeParts[patchi].empty()) continue;

        const polyPatch& patch = pbMesh[patchi];

        const labelList& patchEdgeOwnerOrigBoundaryEdges =
            ncb.patchEdgeOwnerOrigBoundaryEdges(patchi);

        forAll(patchEdgeParts[patchi], patchEdgei)
        {
            const label ownerOrigBoundaryEdgei =
                patchEdgeOwnerOrigBoundaryEdges[patchEdgei];

            const label sign =
                edge::compare
                (
                    meshEdge(patch, patchEdgei),
                    ownerOrigBoundaryMeshEdges[ownerOrigBoundaryEdgei]
                );

            const part& pU = patchEdgeParts[patchi][patchEdgei];
            const part p = sign > 0 ? pU : -pU;

            ownerOrigBoundaryEdgeNParts[ownerOrigBoundaryEdgei] ++;
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei] += p;
        }
    }

    // Give every point in the all boundary a global unique index
    const globalIndex globalOwnerOrigBoundaryPoints
    (
        ownerOrigBoundaryPointMeshPoint.size()
    );
    labelList ownerOrigBoundaryPointIndices
    (
        ownerOrigBoundaryPointMeshPoint.size()
    );
    forAll(ownerOrigBoundaryPointIndices, ownerOrigBoundaryPointi)
    {
        ownerOrigBoundaryPointIndices[ownerOrigBoundaryPointi] =
            globalOwnerOrigBoundaryPoints.toGlobal(ownerOrigBoundaryPointi);
    }
    syncTools::syncPointList
    (
        mesh_,
        ownerOrigBoundaryPointMeshPoint,
        ownerOrigBoundaryPointIndices,
        minEqOp<label>(),
        labelMax
    );

    // Synchronise the sign of the edge parts so that they can be
    // meaningfully summed. This is done using the sign of the global point
    // indices; if they are not in order, then reverse the part area.
    forAll(ownerOrigBoundaryEdgeParts, ownerOrigBoundaryEdgei)
    {
        const edge& e = ownerOrigBoundaryEdges[ownerOrigBoundaryEdgei];

        if
        (
            ownerOrigBoundaryPointIndices[e.start()]
          > ownerOrigBoundaryPointIndices[e.end()])
        {
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei] =
                - ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];
        }
    }

    // Average the boundary edge parts across couplings
    syncTools::syncEdgeList
    (
        mesh_,
        ownerOrigBoundaryEdgeMeshEdge,
        ownerOrigBoundaryEdgeNParts,
        plusEqOp<label>(),
        label(0)
    );
    syncTools::syncEdgeList
    (
        mesh_,
        ownerOrigBoundaryEdgeMeshEdge,
        ownerOrigBoundaryEdgeParts,
        plusEqOp<part>(),
        part(),
        []
        (
            const transformer& vt,
            const bool forward,
            List<part>& fld
        )
        {
            if (forward)
            {
                forAll(fld, i)
                {
                    if (magSqr(fld[i].area) != 0)
                    {
                        fld[i].area = vt.transform(fld[i].area);
                        fld[i].centre = vt.transformPosition(fld[i].centre);
                    }
                }
            }
            else
            {
                forAll(fld, i)
                {
                    if (magSqr(fld[i].area) != 0)
                    {
                        fld[i].area = vt.invTransform(fld[i].area);
                        fld[i].centre = vt.invTransformPosition(fld[i].centre);
                    }
                }
            }
        }
    );
    forAll(ownerOrigBoundaryEdgeParts, ownerOrigBoundaryEdgei)
    {
        if (ownerOrigBoundaryEdgeNParts[ownerOrigBoundaryEdgei] != 0)
        {
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei].area /=
                ownerOrigBoundaryEdgeNParts[ownerOrigBoundaryEdgei];
        }
    }

    // Put the edge parts back to their original orientations
    forAll(ownerOrigBoundaryEdgeParts, ownerOrigBoundaryEdgei)
    {
        const edge& e = ownerOrigBoundaryEdges[ownerOrigBoundaryEdgei];

        if
        (
            ownerOrigBoundaryPointIndices[e.start()]
          > ownerOrigBoundaryPointIndices[e.end()]
        )
        {
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei] =
                - ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];
        }
    }

    return ownerOrigBoundaryEdgeParts;
}


void Foam::fvMeshStitcher::applyOwnerOrigBoundaryEdgeParts
(
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf,
    const List<part>& ownerOrigBoundaryEdgeParts
) const
{
    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList ownerOrigPatchIndices = ncb.ownerOrigPatchIndices();
    const labelList& ownerOrigBoundaryEdgeMeshEdge =
        ncb.ownerOrigBoundaryEdgeMeshEdge();
    const edgeList& ownerOrigBoundaryMeshEdges =
        ncb.ownerOrigBoundaryMeshEdges();

    boolList patchIsOwnerOrig(pbMesh.size(), false);
    UIndirectList<bool>(patchIsOwnerOrig, ownerOrigPatchIndices) = true;

    // Count the number of original faces attached to each boundary edge
    labelList ownerOrigBoundaryEdgeNOrigFaces
    (
        ownerOrigBoundaryEdgeMeshEdge.size(),
        0
    );
    forAll(ownerOrigBoundaryEdgeMeshEdge, ownerOrigBoundaryEdgei)
    {
        const label meshEdgei =
            ownerOrigBoundaryEdgeMeshEdge[ownerOrigBoundaryEdgei];

        forAll(mesh_.edgeFaces()[meshEdgei], edgeFacei)
        {
            const label facei = mesh_.edgeFaces()[meshEdgei][edgeFacei];

            const label patchi =
                mesh_.isInternalFace(facei)
              ? -1
              : pbMesh.patchIndices()[facei - mesh_.nInternalFaces()];

            if (patchi != -1 && patchIsOwnerOrig[patchi])
            {
                ownerOrigBoundaryEdgeNOrigFaces[ownerOrigBoundaryEdgei] ++;
            }
        }
    }

    // Synchronise the boundary edge original face count
    syncTools::syncEdgeList
    (
        mesh_,
        ownerOrigBoundaryEdgeMeshEdge,
        ownerOrigBoundaryEdgeNOrigFaces,
        plusEqOp<label>(),
        label(0)
    );

    // If debugging, check that face changes are in sync
    tmp<surfaceLabelField> tChanged =
        debug
      ? surfaceLabelField::New("changed", mesh_, dimensioned<label>(dimless, 0))
      : tmp<surfaceLabelField>(nullptr);

    // Modify faces connected to the boundary edges
    forAll(ownerOrigBoundaryEdgeMeshEdge, ownerOrigBoundaryEdgei)
    {
        const label meshEdgei =
            ownerOrigBoundaryEdgeMeshEdge[ownerOrigBoundaryEdgei];

        const part& ownerOrigBoundaryEdgeP =
            ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];

        switch (ownerOrigBoundaryEdgeNOrigFaces[ownerOrigBoundaryEdgei])
        {
            case 0:
            {
                continue;
            }

            // If this edge is connected to just one original face then add any
            // edge area into that face
            case 1:
            {
                label origFacei = -1, origPatchi = -1, origPatchFacei = -1;
                forAll(mesh_.edgeFaces()[meshEdgei], edgeFacei)
                {
                    const label facei = mesh_.edgeFaces()[meshEdgei][edgeFacei];

                    const label patchi =
                        mesh_.isInternalFace(facei)
                      ? -1
                      : pbMesh.patchIndices()[facei - mesh_.nInternalFaces()];

                    if (patchi != -1 && patchIsOwnerOrig[patchi])
                    {
                        origFacei = facei;
                        origPatchi = patchi;
                        origPatchFacei = facei - pbMesh[patchi].start();
                        break;
                    }
                }

                if (origFacei != -1)
                {
                    vector& area =
                        SfSf.boundaryFieldRef()[origPatchi][origPatchFacei];
                    point& centre =
                        CfSf.boundaryFieldRef()[origPatchi][origPatchFacei];

                    const label sign =
                        mesh_.faces()[origFacei].edgeDirection
                        (
                            ownerOrigBoundaryMeshEdges[ownerOrigBoundaryEdgei]
                        );

                    part p(area, centre);
                    p +=
                        sign > 0
                      ? ownerOrigBoundaryEdgeP
                      : -ownerOrigBoundaryEdgeP;

                    area = p.area;
                    centre = p.centre;

                    if (debug)
                    {
                        tChanged->boundaryFieldRef()
                            [origPatchi][origPatchFacei] = label(1);
                    }
                }

                break;
            }

            // If this edge is connected to two original faces then add any
            // edge area into non-original connected faces
            case 2:
            {
                forAll(mesh_.edgeFaces()[meshEdgei], edgeFacei)
                {
                    const label facei = mesh_.edgeFaces()[meshEdgei][edgeFacei];

                    const label patchi =
                        mesh_.isInternalFace(facei)
                      ? -1
                      : pbMesh.patchIndices()[facei - mesh_.nInternalFaces()];

                    if (patchi != -1 && patchIsOwnerOrig[patchi])
                    {
                        continue;
                    }

                    const label patchFacei =
                        patchi == -1 ? -1 : facei - pbMesh[patchi].start();

                    vector& area =
                        patchi == -1
                      ? SfSf[facei]
                      : SfSf.boundaryFieldRef()[patchi][patchFacei];
                    point& centre =
                        patchi == -1
                      ? CfSf[facei]
                      : CfSf.boundaryFieldRef()[patchi][patchFacei];

                    const label sign =
                        mesh_.faces()[facei].edgeDirection
                        (
                            ownerOrigBoundaryMeshEdges[ownerOrigBoundaryEdgei]
                        );

                    part p(area, centre);
                    p +=
                        sign < 0
                      ? ownerOrigBoundaryEdgeP
                      : -ownerOrigBoundaryEdgeP;

                    area = p.area;
                    centre = p.centre;

                    if (debug && patchi != -1)
                    {
                        tChanged->boundaryFieldRef()
                            [patchi][patchFacei] = label(1);
                    }
                }

                break;
            }

            default:
            {
                FatalErrorInFunction
                    << "Boundary is not manifold"
                    << exit(FatalError);
            }
        }
    }

    // If debugging, check that face changes are in sync
    if (debug)
    {
        const surfaceLabelField::Boundary& changedBf =
            tChanged->boundaryField();

        const surfaceLabelField::Boundary changedNbrBf
        (
            surfaceLabelField::null(),
            changedBf.boundaryNeighbourField()
        );

        bool synchronised = true;

        forAll(changedBf, patchi)
        {
            const fvsPatchLabelField& changedPf = changedBf[patchi];
            const fvsPatchLabelField& changedNbrPf = changedNbrBf[patchi];

            if (!isA<processorFvPatch>(changedPf.patch())) continue;

            forAll(changedPf, patchFacei)
            {
                if (changedPf[patchFacei] != changedNbrPf[patchFacei])
                {
                    const label facei =
                        changedPf.patch().start() + patchFacei;

                    FatalErrorInFunction
                        << "Face not synchronised with centre at "
                        << mesh_.faceCentres()[facei] << endl;

                    synchronised = false;
                }
            }
        }

        if (!synchronised)
        {
            FatalErrorInFunction
                << exit(FatalError);
        }
    }
}


void Foam::fvMeshStitcher::stabiliseOrigPatchFaces
(
    surfaceVectorField::Boundary& SfBf,
    surfaceVectorField::Boundary& CfBf
) const
{
    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList allOrigPatchIndices = ncb.allOrigPatchIndices();

    forAll(allOrigPatchIndices, i)
    {
        const label origPatchi = allOrigPatchIndices[i];
        const polyPatch& origPp = mesh_.boundaryMesh()[origPatchi];

        forAll(origPp, origPatchFacei)
        {
            const vector& a = origPp.faceAreas()[origPatchFacei];
            const point& c = origPp.faceCentres()[origPatchFacei];

            vector& Sf = SfBf[origPatchi][origPatchFacei];
            point& Cf = CfBf[origPatchi][origPatchFacei];

            // Determine the direction in which to stabilise. If the fv-face
            // points in the same direction as the poly-face, then stabilise in
            // the direction of the poly-face. If it is reversed, then
            // construct a tangent to both faces, and stabilise in the average
            // direction to this tangent and the poly-face.
            vector dSfHat;
            if ((Sf & a) >= 0)
            {
                dSfHat = normalised(a);
            }
            else
            {
                dSfHat = (Sf & Sf)*a - (Sf & a)*Sf;
                if ((dSfHat & a) <= 0) dSfHat = perpendicular(a);
                dSfHat = normalised(normalised(dSfHat) + normalised(a));
            }

            part SAndCf(Sf, Cf);

            SAndCf += part(small*mag(a)*dSfHat, c);

            Sf = SAndCf.area;
            Cf = SAndCf.centre;
        }
    }
}


void Foam::fvMeshStitcher::intersect
(
    surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf,
    const boolList& patchCoupleds,
    const bool matchTopology
) const
{
    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh_);
    const labelList ownerOrigPatchIndices = ncb.ownerOrigPatchIndices();
    const edgeList& ownerOrigBoundaryMeshEdges =
        ncb.ownerOrigBoundaryMeshEdges();

    // Alias the boundary geometry fields
    surfaceVectorField::Boundary& SfBf = SfSf.boundaryFieldRef();
    surfaceVectorField::Boundary& CfBf = CfSf.boundaryFieldRef();

    // Create storage for and initialise the edge parts of source patches
    List<List<part>> patchEdgeParts(mesh_.boundary().size());
    forAll(ownerOrigPatchIndices, i)
    {
        const label origPatchi = ownerOrigPatchIndices[i];

        patchEdgeParts[origPatchi].resize
        (
            mesh_.boundaryMesh()[origPatchi].nEdges(),
            part(Zero)
        );
    }

    // If topology is specified, also create a boundary field of the original
    // patch indices on the neighbouring interface, as well as the
    // corresponding areas and centres for stabilisation purposes
    tmp<surfaceLabelField::Boundary> tOrigFacesNbrBf;
    tmp<surfaceVectorField::Boundary> tOrigSfNbrBf;
    tmp<surfaceVectorField::Boundary> tOrigCfNbrBf;
    if (matchTopology)
    {
        getOrigNbrBfs
        (
            polyFacesBf,
            SfBf,
            CfBf,
            tOrigFacesNbrBf,
            tOrigSfNbrBf,
            tOrigCfNbrBf
        );
    }

    // Create intersection geometry for all the coupled non-conformal patches
    forAll(mesh_.boundary(), patchi)
    {
        if (!patchCoupleds[patchi]) continue;

        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (isA<nonConformalCyclicFvPatch>(fvp))
        {
            const nonConformalCyclicFvPatch& nccFvp =
                refCast<const nonConformalCyclicFvPatch>(fvp);

            intersectNonConformalCyclic
            (
                nccFvp,
                polyFacesBf,
                SfBf,
                CfBf,
                tOrigFacesNbrBf,
                tOrigSfNbrBf,
                tOrigCfNbrBf,
                patchEdgeParts[nccFvp.origPatchIndex()]
            );
        }
        else if (isA<nonConformalProcessorCyclicFvPatch>(fvp))
        {}
        else if (isA<nonConformalMappedWallFvPatch>(fvp))
        {
            const nonConformalMappedWallFvPatch& ncmwFvp =
                refCast<const nonConformalMappedWallFvPatch>(fvp);

            intersectNonConformalMappedWall
            (
                ncmwFvp,
                polyFacesBf,
                SfBf,
                CfBf,
                tOrigFacesNbrBf,
                tOrigSfNbrBf,
                tOrigCfNbrBf,
                patchEdgeParts[ncmwFvp.origPatchIndex()]
            );
        }
        else
        {
            FatalErrorInFunction
                << "Coupled patch type not recognised"
                << exit(FatalError);
        }
    }

    // Create small stabilisation geometry for all the any non-coupled
    // non-conformal patches
    forAll(mesh_.boundary(), patchi)
    {
        if (patchCoupleds[patchi]) continue;

        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const polyPatch& origPp =
            refCast<const nonConformalFvPatch>(fvp).origPatch().patch();

        SfBf[patchi] ==
            vectorField
            (
                small*origPp.faceAreas(),
                polyFacesBf[patchi] - origPp.start()
            );
        CfBf[patchi] ==
            vectorField
            (
                origPp.faceCentres(),
                polyFacesBf[patchi] - origPp.start()
            );
    }

    // Construct the boundary edge geometry
    List<part> ownerOrigBoundaryEdgeParts =
        calculateOwnerOrigBoundaryEdgeParts(patchEdgeParts);

    // Add the difference between patch edge parts and all boundary
    // edge parts to the adjacent patch faces. This is an error part.
    forAll(ownerOrigPatchIndices, i)
    {
        const label origPatchi = ownerOrigPatchIndices[i];
        const polyPatch& origPatch = pbMesh[origPatchi];

        const labelList& origPatchEdgeOwnerOrigBoundaryEdges =
            ncb.patchEdgeOwnerOrigBoundaryEdges(origPatchi);

        forAll(patchEdgeParts[origPatchi], origPatchEdgei)
        {
            const label ownerOrigBoundaryEdgei =
                origPatchEdgeOwnerOrigBoundaryEdges[origPatchEdgei];

            const label sign =
                edge::compare
                (
                    meshEdge(origPatch, origPatchEdgei),
                    ownerOrigBoundaryMeshEdges[ownerOrigBoundaryEdgei]
                );

            part errorP =
                patchEdgeParts[origPatchi][origPatchEdgei];
            errorP -=
                sign > 0
              ? ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei]
              : -ownerOrigBoundaryEdgeParts[ownerOrigBoundaryEdgei];

            forAll(origPatch.edgeFaces()[origPatchEdgei], patchEdgeFacei)
            {
                const label patchFacei =
                    origPatch.edgeFaces()[origPatchEdgei][patchEdgeFacei];

                const label sign =
                    origPatch.localFaces()[patchFacei].edgeDirection
                    (
                        origPatch.edges()[origPatchEdgei]
                    );

                part p
                (
                    SfBf[origPatchi][patchFacei],
                    CfBf[origPatchi][patchFacei]
                );
                p += sign > 0 ? errorP : -errorP;

                SfBf[origPatchi][patchFacei] = p.area;
                CfBf[origPatchi][patchFacei] = p.centre;
            }
        }
    }

    // Use the boundary edge geometry to correct all edge-connected faces
    applyOwnerOrigBoundaryEdgeParts(SfSf, CfSf, ownerOrigBoundaryEdgeParts);

    // Stabilise the original patch faces
    stabiliseOrigPatchFaces(SfBf, CfBf);
}


bool Foam::fvMeshStitcher::disconnectThis
(
    const bool changing,
    const bool geometric
)
{
    if (!stitches()) return false;

    // Determine which patches are coupled
    const boolList patchCoupleds =
        geometric
      ? this->patchCoupleds()
      : boolList(mesh_.boundary().size(), false);

    // Map the non-conformal patch field data to the conformal faces in advance
    // of the non-conformal patches being removed
    if (changing)
    {
        preConformSurfaceFields();
        preConformVolFields();
    }

    // Conform the mesh
    if (mesh_.moving())
    {
        surfaceScalarField phi(mesh_.phi());
        for (label i = 1; i < mesh_.phi().nOldTimes(false); ++ i)
        {
            phi.oldTimeRef(i) == mesh_.phi().oldTime(i);
        }
        conformCorrectMeshPhi(phi);
        mesh_.conform(phi);
    }
    else
    {
        mesh_.conform();
    }

    // Resize all the affected patch fields
    resizePatchFields<SurfaceField>();
    resizePatchFields<VolField>();

    // Create null polyTopoChangeMap
    const polyTopoChangeMap map(mesh_);

    meshObjects::topoChange<fvMesh>(mesh_, map);
    meshObjects::topoChange<lduMesh>(mesh_, map);

    const_cast<Time&>(mesh_.time()).functionObjects().topoChange(map);

    return true;
}


bool Foam::fvMeshStitcher::connectThis
(
    const bool changing,
    const bool geometric,
    const bool load
)
{
    if (!stitches()) return false;

    // Create a copy of the conformal poly face addressing
    IOobject polyFacesBfIO(word::null, mesh_.pointsInstance(), mesh_);
    surfaceLabelField::Boundary polyFacesBf
    (
        surfaceLabelField::null(),
        mesh_.polyFacesBf()
    );

    // If starting up then load topology from disk, if it is available
    const bool matchTopology =
        load && loadPolyFacesBf(polyFacesBfIO, polyFacesBf);

    // Determine which patches are coupled
    const boolList patchCoupleds =
        (geometric || !matchTopology)
      ? this->patchCoupleds()
      : boolList(mesh_.boundary().size(), false);

    // Access all the intersections in advance. Makes the log nicer.
    forAll(mesh_.boundary(), patchi)
    {
        if (!patchCoupleds[patchi]) continue;

        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (isA<nonConformalCyclicFvPatch>(fvp))
        {
            const nonConformalCyclicFvPatch& nccFvp =
                refCast<const nonConformalCyclicFvPatch>(fvp);

            if (nccFvp.owner())
            {
                nccFvp.nonConformalCyclicPatch().intersection();
            }
        }
        else if (isA<nonConformalProcessorCyclicFvPatch>(fvp))
        {}
        else if (isA<nonConformalMappedWallFvPatch>(fvp))
        {
            const nonConformalMappedWallFvPatch& ncmwFvp =
                refCast<const nonConformalMappedWallFvPatch>(fvp);

            const nonConformalMappedWallFvPatch& ownerNcmwFvp =
                ncmwFvp.owner() ? ncmwFvp : ncmwFvp.nbrPatch();

            ownerNcmwFvp.nonConformalMappedWallPatch().intersection();
        }
        else
        {
            FatalErrorInFunction
                << "Coupled patch type not recognised"
                << exit(FatalError);
        }
    }

    if (any(patchCoupleds))
    {
        Info<< indent << typeName << ": Connecting" << incrIndent << endl;
    }

    // Create copies of geometry fields to be modified
    surfaceVectorField Sf(mesh_.Sf().cloneUnSliced()());
    surfaceVectorField Cf(mesh_.Cf().cloneUnSliced()());

    // Construct non-conformal geometry
    intersect(polyFacesBf, Sf, Cf, patchCoupleds, matchTopology);

    // If the mesh is moving then create any additional non-conformal geometry
    // necessary to correct the mesh fluxes
    if (mesh_.moving())
    {
        createNonConformalCorrectMeshPhiGeometry(polyFacesBf, Sf, Cf);
    }

    // Make the mesh non-conformal
    if (mesh_.moving())
    {
        surfaceScalarField phi(mesh_.phi());
        for (label i = 1; i <= mesh_.phi().nOldTimes(false); ++ i)
        {
            phi.oldTimeRef(i) == mesh_.phi().oldTime(i);
        }
        unconformCorrectMeshPhi(polyFacesBf, Sf, Cf, phi);
        mesh_.unconform(polyFacesBf, Sf, Cf, phi);
    }
    else
    {
        mesh_.unconform(polyFacesBf, Sf, Cf);
    }

    // Reset the poly faces instance to that of any loaded topology
    if (load)
    {
        mesh_.setPolyFacesBfInstance(polyFacesBfIO.instance());
    }

    // Resize all the affected patch fields
    resizePatchFields<SurfaceField>();
    resizePatchFields<VolField>();

    // Map the non-conformal patch field data back from the conformal faces and
    // into the new non-conformal patches
    if (changing)
    {
        postUnconformSurfaceFields();
        postUnconformVolFields();
    }

    // Prevent hangs caused by processor cyclic patches using mesh geometry
    mesh_.deltaCoeffs();

    if (any(patchCoupleds))
    {
        const volScalarField::Internal o(openness());
        const scalar gMaxO = gMax(o);
        Info<< indent << "Cell min/average/max openness = "
            << gMin(o) << '/' << gAverage(o) << '/' << gMaxO << endl;
        if (gMaxO > rootSmall)
        {
            FatalErrorInFunction
                << "Maximum openness of " << gMaxO
                << " is not tolerable" << exit(FatalError);
        }

        if (mesh_.moving())
        {
            for (label i = 0; i <= mesh_.phi().nOldTimes(false); ++ i)
            {
                const volScalarField::Internal vce(volumeConservationError(i));
                const scalar gMaxVce = gMax(vce);
                Info<< indent << "Cell min/average/max ";
                for (label j = 0; j < i; ++ j) Info<< "old-";
                Info<< (i ? "time " : "") << "volume conservation error = "
                    << gMin(vce) << '/' << gAverage(vce) << '/' << gMaxVce
                    << endl;
                if (gMaxVce > rootSmall)
                {
                    FatalErrorInFunction
                        << "Maximum volume conservation error of " << gMaxVce
                        << " is not tolerable" << exit(FatalError);
                }
            }
        }

        if (mesh_.moving())
        {
            // Create a boundary field of the imbalance between the mesh fluxes
            // on either side of interfaces
            surfaceScalarField::Boundary mfe
            (
                surfaceScalarField::Internal::null(),
                mesh_.phi().boundaryField()
            );
            mfe += mesh_.phi().boundaryField().boundaryNeighbourField();

            // Determine the number of non-processor patches
            label nNonProcPatches = 0;
            forAll(mesh_.boundary(), patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];

                if (!isA<processorFvPatch>(fvp))
                {
                    if (nNonProcPatches != patchi)
                    {
                        FatalErrorInFunction
                            << "Processor patches do not follow non-processor "
                            << "patches" << exit(FatalError);
                    }

                    nNonProcPatches ++;
                }
            }

            // Work out the min, avg and max error for all non-conformal
            // cyclics, including any errors on referred processor cyclics
            scalarList minMfe(nNonProcPatches, vGreat);
            scalarField sumMfe(nNonProcPatches, 0), nSumMfe(nNonProcPatches, 0);
            scalarList maxMfe(nNonProcPatches, -vGreat);
            forAll(mfe, patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];

                const label nccPatchi =
                    isA<nonConformalCyclicFvPatch>(fvp)
                  ? refCast<const nonConformalCyclicFvPatch>(fvp)
                   .index()
                  : isA<nonConformalProcessorCyclicFvPatch>(fvp)
                  ? refCast<const nonConformalProcessorCyclicFvPatch>(fvp)
                   .referPatchIndex()
                  : -1;

                if (nccPatchi != -1)
                {
                    minMfe[nccPatchi] =
                        min(minMfe[nccPatchi], min(mfe[patchi]));
                    sumMfe[nccPatchi] += sum(mfe[patchi]);
                    nSumMfe[nccPatchi] += mfe[patchi].size();
                    maxMfe[nccPatchi] =
                        max(maxMfe[nccPatchi], max(mfe[patchi]));
                }
            }
            reduce(minMfe, ListOp<minOp<scalar>>());
            reduce(sumMfe, ListOp<sumOp<scalar>>());
            reduce(nSumMfe, ListOp<sumOp<scalar>>());
            reduce(maxMfe, ListOp<maxOp<scalar>>());

            // Report
            forAll(minMfe, patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];

                if
                (
                    isA<nonConformalCyclicFvPatch>(fvp)
                 && refCast<const nonConformalCyclicFvPatch>(fvp).owner()
                 && nSumMfe[patchi] > 0
                )
                {
                    Info<< indent << fvp.name()
                        << " min/average/max mesh flux error = "
                        << minMfe[patchi] << '/'
                        << sumMfe[patchi]/(nSumMfe[patchi] + vSmall) << '/'
                        << maxMfe[patchi] << endl;
                }
            }
        }

        const volScalarField::Internal pvf(projectedVolumeFraction());
        const scalar gMaxPvf = gMax(pvf);
        Info<< indent << "Cell min/average/max projected volume fraction = "
            << gMin(pvf) << '/' << gAverage(pvf) << '/' << gMaxPvf << endl;
        if (gMaxPvf > minWarnProjectedVolumeFraction_)
        {
            WarningInFunction
                << "Maximum projected volume fraction " << gMaxPvf << " may "
                << "cause instability." << nl << indent << "Volumetric "
                << "distortion can be minimised by making the side of the "
                << "interface" << nl << indent << "with smaller or more "
                << "finely layered cells the neighbour." << nl << indent
                << "This is done by specifying this side second to "
                << "createNonConformalCouples;" << nl << indent << "either on "
                << "the command line, or in the patches (and regions) entries "
                << "within" << nl << indent << "the "
                << "createNonConformalCouplesDict." << endl;
        }
    }

    if (any(patchCoupleds))
    {
        Info<< decrIndent;
    }

    // Create null polyTopoChangeMap
    const polyTopoChangeMap map(mesh_);

    meshObjects::topoChange<fvMesh>(mesh_, map);
    meshObjects::topoChange<lduMesh>(mesh_, map);

    const_cast<Time&>(mesh_.time()).functionObjects().topoChange(map);

    return true;
}


void Foam::fvMeshStitcher::preConformSurfaceFields()
{
    #define PreConformSurfaceFields(Type, nullArg) \
        preConformSurfaceFields<Type>();
    FOR_ALL_FIELD_TYPES(PreConformSurfaceFields);
    #undef PreConformSurfaceFields
}


void Foam::fvMeshStitcher::preConformVolFields()
{
    #define PreConformVolFields(Type, nullArg) \
        preConformVolFields<Type>();
    FOR_ALL_FIELD_TYPES(PreConformVolFields);
    #undef PreConformVolFields
}


template<>
void Foam::fvMeshStitcher::postUnconformSurfaceFields<Foam::vector>()
{
    if (mesh_.topoChanged())
    {
        UPtrList<surfaceVectorField> Ufs(mesh_.curFields<surfaceVectorField>());

        forAll(Ufs, i)
        {
            surfaceVectorField& Uf = Ufs[i];

            const volVectorField& U = surfaceToVolVelocity(Uf);

            if (isNull(U)) Uf.clearOldTimes();
        }
    }

    UPtrList<surfaceVectorField> fields(mesh_.fields<surfaceVectorField>());

    forAll(fields, i)
    {
        conformedFvsPatchField<vector>::unconform
        (
            fields[i].boundaryFieldRefNoStoreOldTimes()
        );
    }
}


void Foam::fvMeshStitcher::postUnconformSurfaceFields()
{
    #define PostUnconformSurfaceFields(Type, nullArg) \
        postUnconformSurfaceFields<Type>();
    FOR_ALL_FIELD_TYPES(PostUnconformSurfaceFields);
    #undef PostUnconformSurfaceFields
}


void Foam::fvMeshStitcher::postUnconformVolFields()
{
    #define PostUnconformVolFields(Type, nullArg) \
        postUnconformVolFields<Type>();
    FOR_ALL_FIELD_TYPES(PostUnconformVolFields);
    #undef PostUnconformVolFields

    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh());
    const labelList origPatchIndices = ncb.allOrigPatchIndices();
    const labelList errorPatchIndices = ncb.allErrorPatchIndices();

    // Special handling for velocities. Wherever we find a movingWall-type
    // boundary condition on an original patch, override the corresponding
    // error patch condition to movingWallSlipVelocity.
    UPtrList<volVectorField> Us(mesh().fields<volVectorField>());
    forAll(Us, i)
    {
        volVectorField& U = Us[i];

        forAll(errorPatchIndices, i)
        {
            const label errorPatchi = errorPatchIndices[i];
            const label origPatchi = origPatchIndices[i];

            typename volVectorField::Patch& origUp =
                U.boundaryFieldRefNoStoreOldTimes()[origPatchi];

            if
            (
                isA<movingWallVelocityFvPatchVectorField>(origUp)
             || isA<movingWallSlipVelocityFvPatchVectorField>(origUp)
            )
            {
                U.boundaryFieldRefNoStoreOldTimes().set
                (
                    errorPatchi,
                    new movingWallSlipVelocityFvPatchVectorField
                    (
                        mesh().boundary()[errorPatchi],
                        U
                    )
                );
            }
        }
    }

    #define PostUnconformEvaluateVolFields(Type, nullArg) \
        postUnconformEvaluateVolFields<Type>();
    FOR_ALL_FIELD_TYPES(PostUnconformEvaluateVolFields);
    #undef PostUnconformEvaluateVolFields

    // Special handling for velocities. Recompute the surface velocity using an
    // interpolation of the volume velocity.
    UPtrList<surfaceVectorField> Ufs(mesh_.fields<surfaceVectorField>());
    forAll(Ufs, i)
    {
        surfaceVectorField& Uf = Ufs[i];

        const volVectorField& U = surfaceToVolVelocity(Uf);

        if (isNull(U)) continue;

        const surfaceVectorField UfInterpolated(fvc::interpolate(U));

        forAll(Uf.boundaryField(), patchi)
        {
            if (isA<nonConformalFvPatch>(mesh_.boundary()[patchi]))
            {
                Uf.boundaryFieldRefNoStoreOldTimes()[patchi] ==
                    UfInterpolated.boundaryField()[patchi];
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::boolList Foam::fvMeshStitcher::patchCoupleds() const
{
    boolList result(mesh_.boundary().size(), false);

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        // Initially assume all cyclics can be coupled
        if (isA<nonConformalCyclicFvPatch>(fvp))
        {
            result[patchi] = true;
        }

        // If, however, a processorCyclic is found, and this is not a parallel
        // run, then the associated cyclic cannot be coupled. This works in a
        // single loop because the processorCyclic patches are after the
        // (global) cyclic patches.
        if (isA<nonConformalProcessorCyclicFvPatch>(fvp))
        {
            const nonConformalProcessorCyclicFvPatch& ncpcFvp =
                refCast<const nonConformalProcessorCyclicFvPatch>(fvp);

            result[patchi] = Pstream::parRun();
            result[ncpcFvp.referPatchIndex()] = Pstream::parRun();
        }

        // A mapped wall can be coupled if the neighbour mesh is available and
        // it is a parallel run, or the pair of walls is confined to a single
        // process
        if (isA<nonConformalMappedWallFvPatch>(fvp))
        {
            const nonConformalMappedWallFvPatch& ncmwFvp =
                refCast<const nonConformalMappedWallFvPatch>(fvp);

            result[patchi] =
                ncmwFvp.haveNbr()
             && (
                    Pstream::parRun()
                 || patchToPatchTools::singleProcess
                    (
                        ncmwFvp.patch().size(),
                        ncmwFvp.nbrPatch().patch().size()
                    ) != -1
                );
        }
    }

    return result;
}


bool Foam::fvMeshStitcher::geometric() const
{
    bool result = false;

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!isA<nonConformalFvPatch>(fvp)) continue;

        const fvsPatchScalarField& magSfp =
            mesh_.magSf().boundaryField()[patchi];

        const polyPatch& origPp =
            refCast<const nonConformalFvPatch>(fvp).origPatch().patch();

        const scalarField origMagSfp
        (
            origPp.magFaceAreas(),
            mesh_.polyFacesBf()[patchi] - origPp.start()
        );

        if (max(magSfp/origMagSfp) > rootSmall)
        {
            result = true;
        }
    }

    reduce(result, orOp<bool>());

    return returnReduce(result, orOp<bool>());
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::fvMeshStitcher::openness() const
{
    return
        mag(fvc::surfaceIntegrate(mesh_.Sf()))*mesh_.V()
       /max
        (
            mag(fvc::surfaceSum(cmptMag(mesh_.Sf())))(),
            (small*sqr(cbrt(mesh_.V())))()
        )();
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::fvMeshStitcher::volumeConservationError(const label n) const
{
    if (0 > n || n > 1)
    {
        FatalErrorInFunction
            << "Can only compute volume conservation error for this time, or "
            << "the previous time" << exit(FatalError);
    }

    const surfaceScalarField& phi = mesh_.phi().oldTime(n);

    const dimensionedScalar deltaT =
        n == 0 ? mesh_.time().deltaT() : mesh_.time().deltaT0();

    const volScalarField::Internal& V = n == 0 ? mesh_.V() : mesh_.V0();
    const volScalarField::Internal& V0 = n == 0 ? mesh_.V0() : mesh_.V00();

    return fvc::surfaceIntegrate(phi*deltaT)() - (V - V0)/mesh_.V();
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::fvMeshStitcher::projectedVolumeFraction() const
{
    return mag
    (
        fvc::surfaceIntegrate(mesh_.Sf() & mesh_.Cf())
       /mesh_.nSolutionD()
      - 1
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshStitcher::fvMeshStitcher(fvMesh& mesh)
:
    mesh_(mesh),
    ready_(false),
    regionPolyFacesBfIOs_(),
    regionPolyFacesBfs_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshStitcher::~fvMeshStitcher()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshStitcher::stitches() const
{
    // Meshes in a local sub-directory are not stitched
    if (!mesh_.local().empty())
    {
        return false;
    }

    // Other meshes with non-conformal patches are stitched
    forAll(mesh_.boundary(), patchi)
    {
        if (isA<nonConformalFvPatch>(mesh_.boundary()[patchi]))
        {
            return true;
        }
    }

    return false;
}


bool Foam::fvMeshStitcher::disconnect
(
    const bool changing,
    const bool geometric
)
{
    // Don't do anything if we are already disconnected
    if (mesh_.conformal()) return false;

    // Get all the connected region meshes
    MultiRegionUList<fvMesh> regionMeshes(this->regionMeshes());

    // Disconnect them all
    bool result = false;
    forAll(regionMeshes, i)
    {
        if (regionMeshes[i]().stitcher().disconnectThis(changing, geometric))
        {
            result = true;
        }
    }

    return result;
}


bool Foam::fvMeshStitcher::connect
(
    const bool changing,
    const bool geometric,
    const bool load
)
{
    if (!changing)
    {
        return connectThis(changing, geometric, load);
    }

    // Connection requires all regions to be available and at a consistent
    // state. So, we need to wait until they are all ready. So, if some are not
    // ready, then this function just flips the ready flag for this region. If
    // they are all ready, then this function stitches all of them.

    ready_ = true;

    // Get all the connected region meshes
    MultiRegionUList<fvMesh> regionMeshes(this->regionMeshes());

    // If any mesh is not ready, then don't do anything
    forAll(regionMeshes, i)
    {
        if (!regionMeshes[i]().stitcher().ready_)
        {
            return false;
        }
    }

    // All meshes are ready, so connect all of them up
    bool result = false;
    forAll(regionMeshes, i)
    {
        if (regionMeshes[i]().stitcher().connectThis(changing, geometric, load))
        {
            result = true;
        }
    }

    // Set all meshes as not ready, in preparation for the next iteration
    forAll(regionMeshes, i)
    {
        regionMeshes[i]().stitcher().ready_ = false;
    }

    return result;
}


void Foam::fvMeshStitcher::reconnect(const bool geometric) const
{
    if (mesh_.conformal() || geometric == this->geometric())
    {
        return;
    }

    // Create a copy of the non-conformal poly face addressing
    surfaceLabelField::Boundary polyFacesBf
    (
        surfaceLabelField::null(),
        mesh_.polyFacesBf()
    );

    // Undo all non-conformal changes and clear all geometry and topology
    mesh_.conform();

    // Determine which patches are coupled
    const boolList patchCoupleds =
        geometric
      ? this->patchCoupleds()
      : boolList(mesh_.boundary().size(), false);

    // Create copies of geometry fields to be modified
    surfaceVectorField Sf(mesh_.Sf().cloneUnSliced()());
    surfaceVectorField Cf(mesh_.Cf().cloneUnSliced()());

    // Construct non-conformal geometry
    intersect(polyFacesBf, Sf, Cf, patchCoupleds, true);

    // Apply changes to the mesh
    mesh_.unconform(polyFacesBf, Sf, Cf);

    // Prevent hangs caused by processor cyclic patches using mesh geometry
    mesh_.deltaCoeffs();

    // Create null polyTopoChangeMap
    const polyTopoChangeMap map(mesh_);

    meshObjects::topoChange<fvMesh>(mesh_, map);
    meshObjects::topoChange<lduMesh>(mesh_, map);

    const_cast<Time&>(mesh_.time()).functionObjects().topoChange(map);
}


void Foam::fvMeshStitcher::topoChange(const polyTopoChangeMap&)
{}


void Foam::fvMeshStitcher::mapMesh(const polyMeshMap& map)
{}


void Foam::fvMeshStitcher::distribute(const polyDistributionMap&)
{}


// ************************************************************************* //
