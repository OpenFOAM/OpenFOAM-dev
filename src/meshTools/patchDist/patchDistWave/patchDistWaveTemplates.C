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


template<class PatchPointType, class DataType, class DataMethod>
Foam::label Foam::patchDistWave::getCellValues
(
    const polyMesh& mesh,
    FaceCellWave<PatchPointType>& waveInfo,
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


template
<
    class PatchPointType,
    template<class> class PatchField,
    class DataType,
    class DataMethod
>
Foam::label Foam::patchDistWave::getPatchValues
(
    const polyMesh& mesh,
    FaceCellWave<PatchPointType>& waveInfo,
    FieldField<PatchField, DataType>& patchValues,
    DataMethod method,
    const DataType& stabiliseValue
)
{
    const List<PatchPointType>& faceInfo = waveInfo.allFaceInfo();

    label nInvalid = 0;

    forAll(patchValues, patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        forAll(patchValues[patchi], patchFacei)
        {
            const label facei = patch.start() + patchFacei;

            patchValues[patchi][patchFacei] =
                (faceInfo[facei].*method)(waveInfo.data())
              + stabiliseValue;

            nInvalid += !faceInfo[facei].valid(waveInfo.data());
        }
    }

    return nInvalid;
}


template<class PatchPointType>
Foam::label Foam::patchDistWave::wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance,
    const bool correct
)
{
    // Initialise to changedFacesInfo information to face centre on patches
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
    FaceCellWave<PatchPointType> wave
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceInfo,
        cellInfo,
        mesh.globalData().nTotalCells() + 1 // max iterations
    );

    // Copy distance into return field
    label nUnset = 0;
    scalar (PatchPointType::*dist)(int&) const =
        &PatchPointType::template dist<int>;
    nUnset += getCellValues(mesh, wave, cellDistance, dist);

    // Correct patch cells for true distance
    if (correct)
    {
        Map<label> nearestFace(2*changedFacesInfo.size());
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


template<class PatchPointType, template<class> class PatchField>
Foam::label Foam::patchDistWave::wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    const FieldField<PatchField, typename PatchPointType::dataType>&
        initialPatchData,
    scalarField& cellDistance,
    FieldField<PatchField, scalar>& patchDistance,
    Field<typename PatchPointType::dataType>& cellData,
    FieldField<PatchField, typename PatchPointType::dataType>& patchData,
    const bool correct
)
{
    // Initialise to changedFacesInfo information to face centre on patches
    List<PatchPointType> changedFacesInfo;
    labelList changedFaces;
    setChangedFaces
    (
        mesh,
        patchIDs,
        changedFaces,
        changedFacesInfo,
        initialPatchData
    );

    // Do calculate patch distance by 'growing' from faces.
    List<PatchPointType> faceInfo(mesh.nFaces()), cellInfo(mesh.nCells());
    FaceCellWave<PatchPointType> wave
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceInfo,
        cellInfo,
        mesh.globalData().nTotalCells() + 1 // max iterations
    );

    // Copy distance into return field
    label nUnset = 0;
    scalar (PatchPointType::*dist)(int&) const =
        &PatchPointType::template dist<int>;
    nUnset += getCellValues(mesh, wave, cellDistance, dist);
    nUnset += getPatchValues(mesh, wave, patchDistance, dist, small);

    // Copy data into the return field
    const typename PatchPointType::dataType&
        (PatchPointType::*data)(int&) const =
        &PatchPointType::template data<int>;
    getCellValues(mesh, wave, cellData, data);
    getPatchValues(mesh, wave, patchData, data);

    // Correct patch cells for true distance
    if (correct)
    {
        Map<label> nearestFace(2*changedFacesInfo.size());
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

        // Transfer data from nearest face to cell
        const List<PatchPointType>& faceInfo = wave.allFaceInfo();
        const labelList patchCells(nearestFace.toc());
        forAll(patchCells, patchCelli)
        {
            const label celli = patchCells[patchCelli];
            const label facei = nearestFace[celli];
            cellData[celli] = faceInfo[facei].data();
        }
    }

    return nUnset;
}


// ************************************************************************* //
