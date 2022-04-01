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

#include "patchDistWave.H"
#include "patchDistFuncs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class PatchPointType, class ... InitialPatchData>
void Foam::patchDistWave::setChangedFaces
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    labelList& changedFaces,
    List<PatchPointType>& changedFacesInfo,
    const InitialPatchData& ... initialPatchData
)
{
    label nChangedFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nChangedFaces += mesh.boundaryMesh()[iter.key()].size();
    }

    changedFaces.resize(nChangedFaces);
    changedFacesInfo.resize(nChangedFaces);

    label changedFacei = 0;

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        forAll(patch.faceCentres(), patchFacei)
        {
            const label meshFacei = patch.start() + patchFacei;

            changedFaces[changedFacei] = meshFacei;

            changedFacesInfo[changedFacei] =
                PatchPointType
                (
                    patch.faceCentres()[patchFacei],
                    initialPatchData[patchi][patchFacei] ...,
                    scalar(0)
                );

            changedFacei++;
        }
    }
}


template
<
    class PatchPointType,
    class TrackingData,
    class DataType,
    class DataMethod
>
Foam::label Foam::patchDistWave::getCellValues
(
    FaceCellWave<PatchPointType, TrackingData>& waveInfo,
    Field<DataType>& cellValues,
    DataMethod method,
    const DataType& stabiliseValue
)
{
    const List<PatchPointType>& cellInfo = waveInfo.allCellInfo();

    label nInvalid = 0;

    forAll(cellInfo, celli)
    {
        cellValues[celli] =
            (cellInfo[celli].*method)(waveInfo.data())
          + stabiliseValue;

        nInvalid += !cellInfo[celli].valid(waveInfo.data());
    }

    return nInvalid;
}


template<class PatchPointType, class TrackingData>
Foam::label Foam::patchDistWave::wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance,
    const bool correct,
    TrackingData& td
)
{
    // Initialise changedFacesInfo to face centres on patches
    List<PatchPointType> changedFacesInfo;
    labelList changedFaces;
    setChangedFaces
    (
        mesh,
        patchIDs,
        changedFaces,
        changedFacesInfo
    );

    // Do calculate patch distance by 'growing' from faces.
    List<PatchPointType> faceInfo(mesh.nFaces()), cellInfo(mesh.nCells());
    FaceCellWave<PatchPointType, TrackingData> wave
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceInfo,
        cellInfo,
        mesh.globalData().nTotalCells() + 1, // max iterations
        td
    );

    // Copy distance into return field
    const label nUnset =
        getCellValues
        (
            wave,
            cellDistance,
            &PatchPointType::template dist<TrackingData>
        );

    // Correct patch cells for true distance
    if (correct)
    {
        Map<labelPair> nearestFace(2*changedFacesInfo.size());
        patchDistFuncs::correctBoundaryFaceCells
        (
            mesh,
            patchIDs,
            cellDistance,
            nearestFace
        );
        patchDistFuncs::correctBoundaryPointCells
        (
            mesh,
            patchIDs,
            cellDistance,
            nearestFace
        );
    }

    return nUnset;
}


// ************************************************************************* //
