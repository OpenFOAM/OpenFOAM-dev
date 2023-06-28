/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "fvMeshStitchersMoving.H"
#include "FvFaceCellWave.H"
#include "fvMeshSubset.H"
#include "fvmLaplacian.H"
#include "meshPhiCorrectInfo.H"
#include "meshPhiPreCorrectInfo.H"
#include "movingWallVelocityFvPatchVectorField.H"
#include "movingWallSlipVelocityFvPatchVectorField.H"
#include "regionSplit.H"
#include "solutionControl.H"
#include "syncTools.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
inline List<Type> repeat(const UList<Type>& a, const UList<Type>& b)
{
    List<Type> result(a.size() + b.size());
    forAll(a, i)
    {
        result[2*i] = a[i];
        result[2*i + 1] = b[i];
    }
    return result;
}


template<class Type>
inline List<Type> repeat(const UList<Type>& l)
{
    return repeat(l, l);
}


struct layerAndWeight
{
    label layer;
    scalar weight;

    static const layerAndWeight min, max;

    typedef nil cmptType;

    friend bool operator==(const layerAndWeight& a, const layerAndWeight& b)
    {
        return a.layer == b.layer && a.weight == b.weight;
    }

    friend bool operator!=(const layerAndWeight& a, const layerAndWeight& b)
    {
        return !(a == b);
    }

    friend Ostream& operator<<(Ostream& os, const layerAndWeight& l)
    {
        return os << l.layer << token::SPACE << l.weight;
    }

    friend Istream& operator>>(Istream& is, layerAndWeight& l)
    {
        return is >> l.layer >> l.weight;
    }
};


const layerAndWeight layerAndWeight::min({-labelMax, NaN});


const layerAndWeight layerAndWeight::max({labelMax, NaN});


layerAndWeight max(const layerAndWeight& a, const layerAndWeight& b)
{
    return a.layer > b.layer ? a : b;
}


layerAndWeight min(const layerAndWeight& a, const layerAndWeight& b)
{
    return a.layer < b.layer ? a : b;
}


}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshStitchers
{
    defineTypeNameAndDebug(moving, 0);
    addToRunTimeSelectionTable(fvMeshStitcher, moving, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshStitchers::moving::conformCorrectMeshPhi
(
    surfaceScalarField& phi
)
{
    // Add the non-conformal parts of the mesh flux into the original faces
    surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();
    forAll(phiBf, nccPatchi)
    {
        if (isA<nonConformalFvPatch>(phiBf[nccPatchi].patch()))
        {
            const nonConformalFvPatch& ncFvp =
                refCast<const nonConformalFvPatch>(phiBf[nccPatchi].patch());

            const label origPatchi = ncFvp.origPatchID();
            const fvPatch& origFvp = ncFvp.origPatch();

            const labelList nccOrigPatchFace =
                ncFvp.polyFaces() - origFvp.start();

            for (label i = 0; i <= phi.nOldTimes(); ++ i)
            {
                phi.oldTime(i).boundaryFieldRef()[origPatchi] +=
                    fieldRMapSum
                    (
                        phi.oldTime(i).boundaryField()[nccPatchi],
                        origFvp.size(),
                        nccOrigPatchFace
                    );

                phi.oldTime(i).boundaryFieldRef()[nccPatchi] = 0;
            }
        }
    }
}


void Foam::fvMeshStitchers::moving::createNonConformalCorrectMeshPhiGeometry
(
    surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf
)
{
    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh());
    const labelList origPatchIDs = ncb.allOrigPatchIDs();
    const labelList errorPatchIDs = ncb.allErrorPatchIDs();

    forAll(origPatchIDs, i)
    {
        const label origPatchi = origPatchIDs[i];
        const polyPatch& origPp = mesh().boundaryMesh()[origPatchi];

        const label errorPatchi = errorPatchIDs[i];

        polyFacesBf[errorPatchi] =
            repeat((identityMap(origPp.size()) + origPp.start())());

        SfSf.boundaryFieldRef()[errorPatchi] =
            repeat
            (
                (rootVSmall*origPp.faceNormals())(),
                (-rootVSmall*origPp.faceNormals())()
            );
        CfSf.boundaryFieldRef()[errorPatchi] =
            repeat(origPp.faceCentres());
    }
}


Foam::labelHashSet Foam::fvMeshStitchers::moving::ownerCoupledCellSet()
{
    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh());

    // For every boundary face, count how many owner-orig faces it is connected
    // to (will be 1 or 0)
    labelList bFaceNSet(mesh().nFaces() - mesh().nInternalFaces(), 0);
    forAll(mesh().boundary(), nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (!isA<nonConformalCoupledFvPatch>(fvp)) continue;

        const nonConformalCoupledFvPatch& nccFvp =
            refCast<const nonConformalCoupledFvPatch>(fvp);

        if (!nccFvp.owner()) continue;

        forAll(nccFvp, nccPatchFacei)
        {
            bFaceNSet
            [
                mesh().polyFacesBf()[nccPatchi][nccPatchFacei]
              - mesh().nInternalFaces()
            ] = 1;
        }
    }

    // For every boundary edge, count how many owner-orig faces it is connected
    // to (should be 0, 1 or 2)
    labelList ownerOrigBoundaryEdgeNSet
    (
        ncb.ownerOrigBoundaryEdgeMeshEdge().size(),
        0
    );
    forAll(ncb.ownerOrigBoundaryEdgeMeshEdge(), ownerOrigBoundaryEdgei)
    {
        const label meshEdgei =
            ncb.ownerOrigBoundaryEdgeMeshEdge()[ownerOrigBoundaryEdgei];

        forAll(mesh().edgeFaces()[meshEdgei], edgeFacei)
        {
            const label facei = mesh().edgeFaces()[meshEdgei][edgeFacei];

            if (mesh().isInternalFace(facei)) continue;

            const label bFacei = facei - mesh().nInternalFaces();

            if (bFaceNSet[bFacei])
            {
                ownerOrigBoundaryEdgeNSet[ownerOrigBoundaryEdgei] += 1;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh(),
        ncb.ownerOrigBoundaryEdgeMeshEdge(),
        ownerOrigBoundaryEdgeNSet,
        plusEqOp<label>(),
        label(0)
    );

    // Isolate the cells that are edge connected to the owner patches. This
    // will form the sub-mesh in which the correction is computed.
    labelHashSet set;
    forAll(ncb.ownerOrigBoundaryEdgeMeshEdge(), ownerOrigBoundaryEdgei)
    {
        if (ownerOrigBoundaryEdgeNSet[ownerOrigBoundaryEdgei] < 2) continue;

        const label meshEdgei =
            ncb.ownerOrigBoundaryEdgeMeshEdge()[ownerOrigBoundaryEdgei];

        forAll(mesh().edgeFaces()[meshEdgei], edgeFacei)
        {
            const label facei = mesh().edgeFaces()[meshEdgei][edgeFacei];

            set.insert(mesh().faceOwner()[facei]);

            if (facei < mesh().nInternalFaces())
            {
                set.insert(mesh().faceNeighbour()[facei]);
            }
        }
    }

    return set;
}


void Foam::fvMeshStitchers::moving::unconformInternalFaceCorrectMeshPhi
(
    surfaceScalarField& phi
)
{
    const surfaceScalarField::Boundary& magSfBf =
        mesh().magSf().boundaryField();

    surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

    // Step 1: Construct some boundary information

    // For each boundary face mark if it is an owner orig face and sum the total
    // non-conformal coupled area if so
    boolList bFaceIsOwnerOrig(mesh().nFaces() - mesh().nInternalFaces(), false);
    scalarList bFaceNccMagSf(mesh().nFaces() - mesh().nInternalFaces(), Zero);
    forAll(mesh().boundary(), nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (!isA<nonConformalCoupledFvPatch>(fvp)) continue;

        const nonConformalCoupledFvPatch& nccFvp =
            refCast<const nonConformalCoupledFvPatch>(fvp);

        if (!nccFvp.owner()) continue;

        forAll(nccFvp, nccPatchFacei)
        {
            const label bFacei =
                nccFvp.polyFaces()[nccPatchFacei] - mesh().nInternalFaces();

            bFaceIsOwnerOrig[bFacei] = true;
            bFaceNccMagSf[bFacei] += magSfBf[nccPatchi][nccPatchFacei];
        }
    }

    // Step 2: Construct a sub-mesh for all cells that are connected to the
    // owner patches.

    // Isolate the cells that are edge connected to the owner patches. This
    // will form the sub-mesh in which the correction is computed.
    labelHashSet subCellSet = ownerCoupledCellSet();

    // Create a mesh for the cells that are edge connected to coupled faces
    // of the owner patches. This is where the correction will be computed.
    fvMeshSubset subsetter(mesh());
    subsetter.setLargeCellSubset(ownerCoupledCellSet());
    const fvMesh& subMesh = subsetter.subMesh();
    subMesh.deltaCoeffs();

    // Determine the disconnected regions of the sub mesh
    const regionSplit subMeshRegions(subMesh);
    const label subNRegions = subMeshRegions.nRegions();

    // Map from mesh boundary face to sub-mesh region
    labelList bFaceSubRegion(mesh().nFaces() - mesh().nInternalFaces(), -1);
    forAll(subsetter.faceMap(), subFacei)
    {
        const label facei = subsetter.faceMap()[subFacei];

        if (mesh().isInternalFace(facei)) continue;

        const label bFacei = facei - mesh().nInternalFaces();

        bFaceSubRegion[bFacei] =
            subMeshRegions[subMesh.faceOwner()[subFacei]];
    }

    // Get a single reference cell for each region
    labelList subMeshRegionRefCells(subNRegions, -1);
    {
        static const label proci = Pstream::myProcNo();

        labelList subMeshRegionRefProcs(subNRegions, labelMax);
        forAll(subMeshRegions, subCelli)
        {
            subMeshRegionRefProcs[subMeshRegions[subCelli]] = proci;
        }
        reduce(subMeshRegionRefProcs, ListOp<minOp<label>>());

        forAll(subMeshRegions, subCelli)
        {
            if
            (
                subMeshRegionRefProcs[subMeshRegions[subCelli]] == proci
             && subMeshRegionRefCells[subMeshRegions[subCelli]] == -1
            )
            {
                subMeshRegionRefCells[subMeshRegions[subCelli]] = subCelli;
            }
        }
    }

    // Step 3: Synchronise the flux across the interfaces

    // Create a synchronised mesh flux for the coupled patches. Set both sides
    // to the neighbour value so that all the error is on the owner side.
    surfaceScalarField::Boundary syncPhiBf
    (
        surfaceScalarField::Internal::null(),
        synchronisedBoundaryField<scalar>(phiBf, true, 0, 1)
    );

    // Determine the total mesh flux error and area magnitude for each region
    scalarList regionPhiError(subNRegions, scalar(0));
    scalarList regionMagSf(subNRegions, vSmall);
    forAll(phiBf, nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (!isA<nonConformalCoupledFvPatch>(fvp)) continue;

        const nonConformalCoupledFvPatch& nccFvp =
            refCast<const nonConformalCoupledFvPatch>(fvp);

        if (!nccFvp.owner()) continue;

        forAll(nccFvp, nccPatchFacei)
        {
            const label subRegioni =
                bFaceSubRegion
                [
                    mesh().polyFacesBf()[nccPatchi][nccPatchFacei]
                  - mesh().nInternalFaces()
                ];

            regionPhiError[subRegioni] +=
                syncPhiBf[nccPatchi][nccPatchFacei]
              - phiBf[nccPatchi][nccPatchFacei];

            regionMagSf[subRegioni] +=
                magSfBf[nccPatchi][nccPatchFacei];
        }
    }
    reduce(regionPhiError, ListOp<sumOp<scalar>>());
    reduce(regionMagSf, ListOp<sumOp<scalar>>());

    // Synchronise the mesh fluxes, but offset so that the total flux for each
    // region is the same as for the non-synchronised mesh fluxes
    forAll(phiBf, nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (!isA<nonConformalCoupledFvPatch>(fvp)) continue;

        const nonConformalCoupledFvPatch& nccFvp =
            refCast<const nonConformalCoupledFvPatch>(fvp);

        if (!nccFvp.owner()) continue;

        forAll(nccFvp, nccPatchFacei)
        {
            const label subRegioni =
                bFaceSubRegion
                [
                    mesh().polyFacesBf()[nccPatchi][nccPatchFacei]
                  - mesh().nInternalFaces()
                ];

            phiBf[nccPatchi][nccPatchFacei] =
                syncPhiBf[nccPatchi][nccPatchFacei]
              - magSfBf[nccPatchi][nccPatchFacei]
               /regionMagSf[subRegioni]
               *regionPhiError[subRegioni];
        }
    }

    // Step 4: Set up a system on the sub-mesh with which to solve for a flux
    // that corrects the volume conservation error in the cells connected to
    // synchronised faces

    // Map volumes to the sub mesh
    const volScalarField::Internal subV(subsetter.interpolate(mesh().V()));
    const volScalarField::Internal subV0(subsetter.interpolate(mesh().V0()));

    // Map mesh flux to the sub mesh, accumulating when a face is split into
    // multiple non-conformal parts
    surfaceScalarField subPhi
    (
        surfaceScalarField::New
        (
            "phi",
            subMesh,
            dimensionedScalar(dimVolume/dimTime, 0)
        )
    );
    forAll(subPhi, subFacei)
    {
        const label facei = subsetter.faceMap()[subFacei];

        subPhi[subFacei] = phi[facei];
    }
    forAll(subPhi.boundaryField(), subPatchi)
    {
        const fvPatch& subFvp = subPhi.boundaryField()[subPatchi].patch();

        forAll(subPhi.boundaryField()[subPatchi], subPatchFacei)
        {
            const label subFacei = subFvp.start() + subPatchFacei;
            const label facei = subsetter.faceMap()[subFacei];

            if (mesh().isInternalFace(facei))
            {
                const label s = sign(subsetter.faceFlipMap()[subFacei]);

                subPhi.boundaryFieldRef()[subPatchi][subPatchFacei] =
                    s*phi[facei];
            }
            else
            {
                const label bFacei = facei - mesh().nInternalFaces();

                const labelUList patches =
                    mesh().polyBFacePatches()[bFacei];
                const labelUList patchFaces =
                    mesh().polyBFacePatchFaces()[bFacei];

                forAll(patches, i)
                {
                    subPhi.boundaryFieldRef()[subPatchi][subPatchFacei] +=
                        phiBf[patches[i]][patchFaces[i]];
                }
            }
        }
    }

    // Calculate the volume conservation error for the sub mesh
    const volScalarField::Internal subVce
    (
        fvc::surfaceIntegrate(subPhi*subMesh.time().deltaT())()
      - (subV - subV0)/subV
    );

    // Construct boundary conditions for the sub-mesh potential. Zero gradient
    // is used for all non-constraint boundaries (i.e., this is a closed
    // domain) and for exposed internal faces, which therefore require a patch
    // type override
    wordList MeshPhiPatchTypes(subMesh.boundary().size());
    forAll(subMesh.boundary(), patchi)
    {
        const fvPatch& subFvp = subMesh.boundary()[patchi];
        MeshPhiPatchTypes[patchi] =
            polyPatch::constraintType(subFvp.type())
         && !isA<internalFvPatch>(subFvp)
          ? subFvp.type()
          : zeroGradientFvPatchField<scalar>::typeName;
    }

    // Solve
    volScalarField MeshPhi
    (
        IOobject
        (
            "MeshPhi",
            subMesh.time().name(),
            subMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        subMesh,
        dimensionedScalar(dimArea, 0),
        MeshPhiPatchTypes,
        subMesh.boundaryMesh().types()
    );

    subMesh.schemes().setFluxRequired(MeshPhi.name());

    fvScalarMatrix MeshPhiEqn
    (
        fvm::laplacian(MeshPhi) + subVce
    );

    forAll(subMeshRegionRefCells, i)
    {
        MeshPhiEqn.setReference(subMeshRegionRefCells[i], 0);
    }

    MeshPhiEqn.solve();

    // Create a correction for the sub-mesh face fluxes
    surfaceScalarField subDeltaPhi(MeshPhiEqn.flux()/mesh().time().deltaT());

    // Step 5: Currently, the computed correction fixes the volume conservation
    // error to the level of the solve tolerance. We would rather the volume
    // conservation error be at round off and the error be in the mismatch
    // across the couplings, as the mismatch will be fixed later when area and
    // mesh flux is added to the error faces. So, here, we correct the
    // correction (!) by propagating the error from the cells furthest from the
    // interface back to the interface. This is done using waves.

    // Dynamic memory for wave initialisation
    DynamicList<labelPair> subChangedPatchAndFaces;
    DynamicList<meshPhiPreCorrectInfo> subChangedFacePci;
    DynamicList<meshPhiCorrectInfo> subChangedFaceCi;

    // Allocate data for the pre-correction wave
    List<meshPhiPreCorrectInfo> subInternalFacePci(subMesh.nInternalFaces());
    List<List<meshPhiPreCorrectInfo>> subPatchFacePci
    (
        FvFaceCellWave<meshPhiPreCorrectInfo>::template
        sizesListList<List<List<meshPhiPreCorrectInfo>>>
        (
            FvFaceCellWave<meshPhiPreCorrectInfo>::template
            listListSizes(subMesh.boundary()),
            meshPhiPreCorrectInfo()
        )
    );
    List<meshPhiPreCorrectInfo> subCellPci(subMesh.nCells());

    // Initialisation for the pre-correction wave
    subChangedPatchAndFaces.clear();
    subChangedFacePci.clear();
    forAll(subMesh.boundary(), subPatchi)
    {
        const fvPatch& subFvp = subMesh.boundary()[subPatchi];

        forAll(subFvp, subPatchFacei)
        {
            const label subFacei = subFvp.start() + subPatchFacei;
            const label facei = subsetter.faceMap()[subFacei];
            const label bFacei = facei - mesh().nInternalFaces();

            if (bFacei >= 0 && bFaceIsOwnerOrig[bFacei])
            {
                subChangedPatchAndFaces.append({subPatchi, subPatchFacei});
                subChangedFacePci.append
                (
                    meshPhiPreCorrectInfo(0, bFaceNccMagSf[bFacei])
                );
            }
        }
    }

    // Pre-correction wave
    FvFaceCellWave<meshPhiPreCorrectInfo> preWave
    (
        subMesh,
        subInternalFacePci,
        subPatchFacePci,
        subCellPci
    );
    preWave.setFaceInfo(subChangedPatchAndFaces, subChangedFacePci);
    const label nWaveLayers =
        preWave.iterate(subMesh.globalData().nTotalCells() + 1);

    // Allocate data for the correction wave
    List<meshPhiCorrectInfo> subInternalFaceCi(subMesh.nInternalFaces());
    List<List<meshPhiCorrectInfo>> subPatchFaceCi
    (
        FvFaceCellWave<meshPhiCorrectInfo>::template
        sizesListList<List<List<meshPhiCorrectInfo>>>
        (
            FvFaceCellWave<meshPhiCorrectInfo>::template
            listListSizes(subMesh.boundary()),
            meshPhiCorrectInfo()
        )
    );
    List<meshPhiCorrectInfo> subCellCi(subMesh.nCells());

    // Calculate the current error in the rate of change of volume
    const volScalarField::Internal dVdtError
    (
        (subV - subV0)/subMesh.time().deltaT()
      - fvc::surfaceIntegrate(subPhi + subDeltaPhi)()*subV
    );

    // Construct track data for the correction wave
    meshPhiCorrectInfo::trackData td
    (
        subInternalFacePci,
        subPatchFacePci,
        subCellPci,
        dVdtError
    );

    // Wave backwards through the layers to generate the corrections. Note that
    // this has to be done in stages, so that later layers complete in their
    // entirety before the earlier layers begin. Otherwise they interfere.
    for (label waveLayeri = nWaveLayers - 1; waveLayeri >= 0; waveLayeri --)
    {
        // The layer indices on the faces that we want to wave from
        const label faceLayeri = (waveLayeri + 1)*2;

        // Initialisation for the correction wave
        subChangedPatchAndFaces.clear();
        subChangedFaceCi.clear();
        forAll(subInternalFacePci, subFacei)
        {
            if (subInternalFacePci[subFacei].layer() == faceLayeri)
            {
                subChangedPatchAndFaces.append({-1, subFacei});
                subChangedFaceCi.append
                (
                    subInternalFaceCi[subFacei].valid(td)
                  ? subInternalFaceCi[subFacei]
                  : meshPhiCorrectInfo(meshPhiCorrectInfo::shape::face)
                );
            }
        }
        forAll(subPatchFacePci, subPatchi)
        {
            forAll(subPatchFacePci[subPatchi], subPatchFacei)
            {
                if
                (
                    subPatchFacePci[subPatchi][subPatchFacei].layer()
                 == faceLayeri
                )
                {
                    subChangedPatchAndFaces.append({subPatchi, subPatchFacei});
                    subChangedFaceCi.append
                    (
                        subPatchFaceCi[subPatchi][subPatchFacei].valid(td)
                      ? subPatchFaceCi[subPatchi][subPatchFacei]
                      : meshPhiCorrectInfo(meshPhiCorrectInfo::shape::face)
                    );
                }
            }
        }

        // Correction wave
        FvFaceCellWave<meshPhiCorrectInfo, meshPhiCorrectInfo::trackData> wave
        (
            subMesh,
            subInternalFaceCi,
            subPatchFaceCi,
            subCellCi,
            td
        );
        wave.setFaceInfo(subChangedPatchAndFaces, subChangedFaceCi);
        wave.iterate(1);
    }

    // Apply corrections
    forAll(subInternalFaceCi, subFacei)
    {
        subDeltaPhi.primitiveFieldRef()[subFacei] +=
            subInternalFaceCi[subFacei].deltaPhi();
    }
    forAll(subPatchFacePci, subPatchi)
    {
        forAll(subPatchFacePci[subPatchi], subPatchFacei)
        {
            subDeltaPhi.boundaryFieldRef()[subPatchi][subPatchFacei] +=
                subPatchFaceCi[subPatchi][subPatchFacei].deltaPhi();
        }
    }

    // Step 6: Apply the corrections to the mesh flux

    // Correct the internal mesh face fluxes
    forAll(subPhi, subFacei)
    {
        phi[subsetter.faceMap()[subFacei]] =
            subPhi[subFacei] + subDeltaPhi[subFacei];
    }

    // Map the sub-mesh flux changes to the conformal mesh boundary
    surfaceScalarField::Boundary deltaPhiBf
    (
        mesh().boundary(),
        surfaceScalarField::Internal::null(),
        calculatedFvPatchField<scalar>::typeName
    );
    deltaPhiBf = 0;
    forAll(subMesh.boundary(), subPatchi)
    {
        const label patchi = subsetter.patchMap()[subPatchi];

        if (patchi == -1) continue;

        const fvPatch& subFvp = subMesh.boundary()[subPatchi];
        const fvPatch& fvp = mesh().boundary()[patchi];

        const bool coupled = subFvp.coupled();

        forAll(subMesh.boundary()[subPatchi], subPatchFacei)
        {
            const label facei =
                subsetter.faceMap()[subFvp.start() + subPatchFacei];

            const label patchFacei = facei - fvp.start();

            if (coupled)
            {
                deltaPhiBf[patchi][patchFacei] =
                    subPhi.boundaryField()[subPatchi][subPatchFacei]
                  + subDeltaPhi.boundaryField()[subPatchi][subPatchFacei]
                  - phiBf[patchi][patchFacei];
            }
            else
            {
                deltaPhiBf[patchi][patchFacei] =
                    subDeltaPhi.boundaryField()[subPatchi][subPatchFacei];
            }
        }
    }

    // Move the changes from the owner-orig to the non-conformal coupled faces
    forAll(mesh().boundary(), nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (!isA<nonConformalCoupledFvPatch>(fvp)) continue;

        const nonConformalCoupledFvPatch& nccFvp =
            refCast<const nonConformalCoupledFvPatch>(fvp);

        if (!nccFvp.owner()) continue;

        const label origPatchi = nccFvp.origPatchID();
        const fvPatch& origFvp = nccFvp.origPatch();

        forAll(nccFvp, nccPatchFacei)
        {
            const label bFacei =
                nccFvp.polyFaces()[nccPatchFacei] - mesh().nInternalFaces();

            const label origPatchFacei =
                nccFvp.polyFaces()[nccPatchFacei] - origFvp.start();

            const scalar deltaPhi =
                magSfBf[nccPatchi][nccPatchFacei]
               /bFaceNccMagSf[bFacei]
               *deltaPhiBf[origPatchi][origPatchFacei];

            deltaPhiBf[nccPatchi][nccPatchFacei] = deltaPhi;
        }
    }
    forAll(mesh().boundary(), nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (!isA<nonConformalCoupledFvPatch>(fvp)) continue;

        const nonConformalCoupledFvPatch& nccFvp =
            refCast<const nonConformalCoupledFvPatch>(fvp);

        if (!nccFvp.owner()) continue;

        const label origPatchi = nccFvp.origPatchID();
        const fvPatch& origFvp = nccFvp.origPatch();

        forAll(nccFvp, nccPatchFacei)
        {
            const label origPatchFacei =
                nccFvp.polyFaces()[nccPatchFacei] - origFvp.start();

            deltaPhiBf[origPatchi][origPatchFacei] = 0;
        }
    }

    // Correct the boundary mesh face fluxes
    phiBf += deltaPhiBf;
}


void Foam::fvMeshStitchers::moving::unconformErrorFaceCorrectMeshPhi
(
    const surfaceLabelField::Boundary& polyFacesBf,
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf,
    surfaceScalarField& phi
)
{
    const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh());
    const labelList origPatchIDs = ncb.allOrigPatchIDs();
    const labelList errorPatchIDs = ncb.allErrorPatchIDs();

    // Synchronise the mesh fluxes on both sides of coupled patches. Store
    // the change made to the mesh flux as an error.
    PtrList<surfaceScalarField::Boundary> phiErrorbs(phi.nOldTimes() + 1);
    for (label i = 0; i <= phi.nOldTimes(); ++ i)
    {
        tmp<surfaceScalarField::Boundary> tphib =
            synchronisedBoundaryField<scalar>(phi.oldTime(i).boundaryField());

        phiErrorbs.set
        (
            i,
            new surfaceScalarField::Boundary
            (
                surfaceScalarField::Internal::null(),
                phi.oldTime(i).boundaryField()
            )
        );
        phiErrorbs[i] = phi.oldTime(i).boundaryField() - tphib();

        phi.oldTime(i).boundaryFieldRef() = tphib;
    }

    // Add the mesh flux error into the error patch so that the mesh fluxes
    // and volume changes still match
    forAll(mesh().boundary(), nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (isA<nonConformalCoupledFvPatch>(fvp))
        {
            const nonConformalCoupledFvPatch& nccFvp =
                refCast<const nonConformalCoupledFvPatch>(fvp);

            const label origPatchi = nccFvp.origPatchID();
            const polyPatch& origPp = mesh().boundaryMesh()[origPatchi];

            const label errorPatchi = nccFvp.errorPatchID();

            forAll(nccFvp, nccPatchFacei)
            {
                const label origPatchFacei =
                    nccFvp.polyFaces()[nccPatchFacei] - origPp.start();

                const label errorPatchFacei0 = 2*origPatchFacei;
                const label errorPatchFacei1 = 2*origPatchFacei + 1;

                for (label i = 0; i <= phi.nOldTimes(); ++ i)
                {
                    fvsPatchField<scalar>& phip =
                        phi.oldTime(i).boundaryFieldRef()[errorPatchi];
                    phip[errorPatchFacei0] +=
                        phiErrorbs[i][nccPatchi][nccPatchFacei]/2;
                    phip[errorPatchFacei1] +=
                        phiErrorbs[i][nccPatchi][nccPatchFacei]/2;
                }
            }
        }
    }

    // Create a boundary field with a mesh velocity magnitude. Take the
    // maximum mesh velocity on either side of non-conformal-coupled faces,
    // and then average that into the original faces. This means even
    // stationary original faces have a velocity magnitude stored that is
    // representative of the interface motion.
    tmp<surfaceScalarField> tnccMeshMagUf =
        surfaceScalarField::New
        (
            "nccMeshMagUf",
            mesh(),
            dimensionedScalar(dimVelocity, Zero)
        );
    surfaceScalarField::Boundary& tnccMeshMagUfb =
        tnccMeshMagUf.ref().boundaryFieldRef();
    forAll(mesh().boundary(), nccPatchi)
    {
        const fvPatch& fvp = mesh().boundary()[nccPatchi];

        if (isA<nonConformalCoupledFvPatch>(fvp))
        {
            const nonConformalCoupledFvPatch& nccFvp =
                refCast<const nonConformalCoupledFvPatch>(fvp);

            const fvPatch& origFvp = nccFvp.origPatch();

            forAll(nccFvp, nccPatchFacei)
            {
                const label origPatchFacei =
                    nccFvp.polyFaces()[nccPatchFacei]
                  - origFvp.start();

                const point& origC =
                    origFvp.patch().faceCentres()[origPatchFacei];
                const point origC0 =
                    origFvp.patch()[origPatchFacei]
                   .centre(mesh().oldPoints());

                tnccMeshMagUfb[nccPatchi][nccPatchFacei] =
                    mag(origC - origC0)/mesh().time().deltaTValue();
            }
        }
    }
    tnccMeshMagUfb =
        max
        (
            tnccMeshMagUfb,
            tnccMeshMagUfb.boundaryNeighbourField()()
        );
    surfaceScalarField::Boundary meshMagUfb
    (
        surfaceScalarField::Internal::null(),
        conformalNccBoundaryField<scalar>(tnccMeshMagUfb)
    );
    tnccMeshMagUf.clear();

    // Resize the error patch faces so that mesh flux divided by area
    // results in a velocity equal to the mesh velocity of the original
    // patch face
    forAll(origPatchIDs, i)
    {
        const label origPatchi = origPatchIDs[i];
        const polyPatch& origPp = mesh().boundaryMesh()[origPatchi];

        const label errorPatchi = errorPatchIDs[i];

        forAll(origPp, origPatchFacei)
        {
            const label errorPatchFacei0 = 2*origPatchFacei;
            const label errorPatchFacei1 = 2*origPatchFacei + 1;

            const vector errorSf =
                min
                (
                    mag(phi.boundaryField()[errorPatchi][errorPatchFacei0])
                   /max(meshMagUfb[origPatchi][origPatchFacei], vSmall),
                    origPp.magFaceAreas()[origPatchFacei]
                )
               *origPp.faceNormals()[origPatchFacei];

            fvsPatchField<vector>& Sfp =
                SfSf.boundaryFieldRef()[errorPatchi];
            Sfp[errorPatchFacei0] += errorSf;
            Sfp[errorPatchFacei1] -= errorSf;
        }
    }

    // Wherever we find a movingWall-type boundary condition on an original
    // patch, override the corresponding error patch condition to
    // movingWallSlipVelocity
    UPtrList<volVectorField> fields(mesh().fields<volVectorField>());
    forAll(fields, i)
    {
        volVectorField& field = fields[i];

        typename volVectorField::Boundary& Ub = field.boundaryFieldRef();

        forAll(origPatchIDs, i)
        {
            const label origPatchi = origPatchIDs[i];

            typename volVectorField::Patch& origUp = Ub[origPatchi];

            if
            (
                isA<movingWallVelocityFvPatchVectorField>(origUp)
             || isA<movingWallSlipVelocityFvPatchVectorField>(origUp)
            )
            {
                const label errorPatchi = errorPatchIDs[i];

                Ub.set
                (
                    errorPatchi,
                    new movingWallSlipVelocityFvPatchVectorField
                    (
                        mesh().boundary()[errorPatchi],
                        field
                    )
                );
            }
        }
    }
}


void Foam::fvMeshStitchers::moving::unconformCorrectMeshPhi
(
    const SurfaceFieldBoundary<label>& polyFacesBf,
    surfaceVectorField& SfSf,
    surfaceVectorField& CfSf,
    surfaceScalarField& phi
)
{
    // !!! At present, the correction procedures that follow require the mesh
    // to be unconformed to its final topological state. This means we have to
    // call fvMesh::unconform twice; once to change the topology to allow for
    // calculation of the error areas and mesh flux corrections, and once again
    // to actually apply these error areas and mesh flux corrections. This is
    // OK for now (the operation is relatively quick), but ideally everything
    // would just be calculated in advance and applied to the mesh in one go,
    // and this function would only modify its arguments and leave calling
    // fvMesh::unconform to the base class.
    mesh().unconform(polyFacesBf, SfSf, CfSf);
    resizeFieldPatchFields(polyFacesBf, phi);

    // Set mesh fluxes on the original and cyclic faces as a proportion of
    // the area taken from the old original faces
    for (label i = 0; i <= phi.nOldTimes(); ++ i)
    {
        phi.oldTime(i).boundaryFieldRef() =
            nonConformalBoundaryField<scalar>
            (
                phi.oldTime(i).boundaryField(),
                phi.oldTime(i).boundaryField()
            );
    }

    // Correct the mesh flux error by modifying values on the internal
    // faces that are edge-connected to the owner orig patches.
    //
    // !!! This only corrects the new-time flux. For smooth operation with
    // second order time schemes on wonky meshes, this will probably need to
    // correct old time mesh fluxes, too.
    //
    if
    (
        mesh().foundObject<solutionControl>(solutionControl::typeName)
     && mesh().lookupObject<solutionControl>(solutionControl::typeName)
       .dict().lookup<Switch>("correctMeshPhi")
    )
    {
        unconformInternalFaceCorrectMeshPhi(phi);
    }

    // Correct the mesh flux error by adding flux and area to the error
    // patch faces
    unconformErrorFaceCorrectMeshPhi(polyFacesBf, SfSf, CfSf, phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshStitchers::moving::moving(fvMesh& mesh)
:
    fvMeshStitcher(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshStitchers::moving::~moving()
{}


// ************************************************************************* //
