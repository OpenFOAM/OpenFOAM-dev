/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "fvRegionSplit.H"
#include "FvFaceCellWave.H"
#include "globalIndex.H"
#include "minData.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvRegionSplit, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvRegionSplit::calcNonCompactRegionSplit
(
    const globalIndex& globalFaces,
    const volScalarField& isoField,
    const scalar isoValue,
    labelList& cellRegion
) const
{
    // Initialise wave data. Blocked faces are given a negative number, so it
    // is always less than that of the global face indices, so the walk will
    // not cross them
    List<minData> internalFaceData(mesh().nInternalFaces());
    label nUnblocked = 0;
    if (isNull(isoField))
    {
        nUnblocked = mesh().nInternalFaces();
    }
    else
    {
        forAll(internalFaceData, internalFacei)
        {
            const label own = mesh().owner()[internalFacei];
            const label nbr = mesh().neighbour()[internalFacei];

            if ((isoField[own] < isoValue) != (isoField[nbr] < isoValue))
            {
                internalFaceData[internalFacei] = minData(-2);
            }
            else
            {
                nUnblocked ++;
            }
        }
    }

    List<List<minData>> patchFaceData(mesh().boundary().size());
    forAll(patchFaceData, patchi)
    {
        const fvPatchScalarField& pf = isoField.boundaryField()[patchi];

        patchFaceData[patchi].resize(pf.size());

        if (isNull(isoField) || !pf.coupled())
        {
            nUnblocked += pf.size();
        }
        else
        {
            tmp<scalarField> tpnf =
                refCast<const coupledFvPatchScalarField>(pf)
               .patchNeighbourField();
            const scalarField& pnf = tpnf();

            forAll(patchFaceData[patchi], patchFacei)
            {
                if ((pf[patchFacei] < isoValue) != (pnf[patchFacei] < isoValue))
                {
                    patchFaceData[patchi][patchFacei] = minData(-2);
                }
                else
                {
                    nUnblocked ++;
                }
            }
        }
    }

    List<minData> cellData(mesh().nCells());

    // Seed all unblocked faces with a globally unique number
    List<labelPair> seedPatchAndFaces(nUnblocked);
    List<minData> seedData(nUnblocked);
    auto& td = FvFaceCellWave<minData>::defaultTrackingData_;
    nUnblocked = 0;
    forAll(internalFaceData, internalFacei)
    {
        if (!internalFaceData[internalFacei].valid(td))
        {
            seedPatchAndFaces[nUnblocked] = labelPair(-1, internalFacei);
            seedData[nUnblocked] = minData(globalFaces.toGlobal(internalFacei));
            nUnblocked ++;
        }
    }
    forAll(patchFaceData, patchi)
    {
        forAll(patchFaceData[patchi], patchFacei)
        {
            if (!patchFaceData[patchi][patchFacei].valid(td))
            {
                seedPatchAndFaces[nUnblocked] = labelPair(patchi, patchFacei);
                seedData[nUnblocked] =
                    minData
                    (
                        globalFaces.toGlobal
                        (
                            mesh().polyFacesBf()[patchi][patchFacei]
                        )
                    );
                nUnblocked ++;
            }
        }
    }

    // Propagate information inwards
    FvFaceCellWave<minData> deltaCalc
    (
        mesh(),
        seedPatchAndFaces,
        seedData,
        internalFaceData,
        patchFaceData,
        cellData,
        mesh().globalData().nTotalCells() + 1,
        FvFaceCellWave<minData>::defaultTrackingData_
    );

    // Extract into the cells
    cellRegion.setSize(mesh().nCells());
    forAll(cellRegion, celli)
    {
        if (cellData[celli].valid(deltaCalc.data()))
        {
            cellRegion[celli] = cellData[celli].data();
        }
        else
        {
            // Unvisited cell. This is only possible if it is surrounded by
            // blocked faces. If so make up region from any of the faces.
            const label facei = mesh().cells()[celli][0];
            cellRegion[celli] = globalFaces.toGlobal(facei);
        }
    }
}


Foam::label Foam::fvRegionSplit::calcRegionSplit
(
    const volScalarField& isoField,
    const scalar isoValue,
    labelList& cellRegion
) const
{
    // Create global face index
    const globalIndex globalFaces(mesh().nFaces());

    // Minimise regions across connected cells. This leaves cellRegion
    // containing the minimum face index in the contiguous region of mesh.
    calcNonCompactRegionSplit(globalFaces, isoField, isoValue, cellRegion);

    // Compact and return
    return compactGlobalRegionSplit(globalFaces, cellRegion);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvRegionSplit::fvRegionSplit(const fvMesh& mesh)
:
    regionSplitBase(mesh.nCells()),
    mesh_(mesh)
{
    nRegions_ =
        calcRegionSplit
        (
            NullObjectRef<volScalarField>(),
            NaN,
            *this
        );
}


Foam::fvRegionSplit::fvRegionSplit
(
    const fvMesh& mesh,
    const volScalarField& isoField,
    const scalar isoValue
)
:
    regionSplitBase(mesh.nCells()),
    mesh_(mesh)
{
    nRegions_ =
        calcRegionSplit
        (
            isoField,
            isoValue,
            *this
        );
}


// ************************************************************************* //
