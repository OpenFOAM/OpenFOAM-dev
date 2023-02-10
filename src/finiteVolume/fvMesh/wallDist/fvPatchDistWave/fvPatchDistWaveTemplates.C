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

#include "fvPatchDistWave.H"
#include "FvWallInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvPatchDistWave
{

template<class FvWallInfoType, class TrackingData>
const List<FvWallInfoType>& getInternalInfo
(
    const volScalarField& distance,
    FvFaceCellWave<FvWallInfoType, TrackingData>& wave
)
{
    return wave.cellInfo();
}

template<class FvWallInfoType, class TrackingData>
const List<FvWallInfoType>& getInternalInfo
(
    const surfaceScalarField& distance,
    FvFaceCellWave<FvWallInfoType, TrackingData>& wave
)
{
    return wave.internalFaceInfo();
}

}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class FvWallInfoType,
    class TrackingData,
    template<class> class PatchField,
    class GeoMesh,
    class ... DataType
>
Foam::label Foam::fvPatchDistWave::wave
(
    const fvMesh& mesh,
    const List<labelPair>& changedPatchAndFaces,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    TrackingData& td,
    GeometricField<DataType, PatchField, GeoMesh>& ... data
)
{
    // If the number of corrections is less than 0 (i.e., -1) then this is a
    // calculation across the entire mesh. Otherwise it is a correction for the
    // cells/faces near the changed faces.
    const bool calculate = nCorrections < 0;

    // Quick return if no corrections
    if (!calculate && nCorrections == 0) return 0;

    // Initialise changedFacesInfo to face centres on patches
    List<FvWallInfoType> changedFacesInfo(changedPatchAndFaces.size());
    forAll(changedPatchAndFaces, changedFacei)
    {
        const label patchi =
            changedPatchAndFaces[changedFacei].first();
        const label patchFacei =
            changedPatchAndFaces[changedFacei].second();

        changedFacesInfo[changedFacei] =
            FvWallInfoType
            (
                data.boundaryField()[patchi][patchFacei] ...,
                mesh.boundaryMesh()[patchi][patchFacei],
                mesh.points(),
                mesh.Cf().boundaryField()[patchi][patchFacei],
                scalar(0)
            );
    }

    // Do calculate patch distance by 'growing' from faces.
    List<FvWallInfoType> internalFaceInfo(mesh.nInternalFaces());
    List<List<FvWallInfoType>> patchFaceInfo
    (
        FvFaceCellWave<FvWallInfoType, TrackingData>::template
        sizesListList<List<List<FvWallInfoType>>>
        (
            FvFaceCellWave<FvWallInfoType, TrackingData>::template
            listListSizes(mesh.boundary()),
            FvWallInfoType()
        )
    );
    List<FvWallInfoType> cellInfo(mesh.nCells());

    // Prevent hangs associated with generation of on-demand geometry
    mesh.C();
    mesh.Cf();

    // Do the wave
    FvFaceCellWave<FvWallInfoType, TrackingData> wave
    (
        mesh,
        internalFaceInfo,
        patchFaceInfo,
        cellInfo,
        td
    );
    wave.setFaceInfo(changedPatchAndFaces, changedFacesInfo);
    if (calculate)
    {
        // Calculation. Wave to completion.
        wave.iterate(mesh.globalData().nTotalCells() + 1);
    }
    else
    {
        // Correction. Wave the specified number of times then stop. We care
        // about cell values, so avoid the final cellToFace by doing n - 1
        // iterations than a final faceToCell.
        wave.iterate(nCorrections - 1);
        wave.faceToCell();
    }

    // Copy distances into field
    const List<FvWallInfoType>& internalInfo = getInternalInfo(distance, wave);
    label nUnset = 0;
    forAll(internalInfo, internali)
    {
        const bool valid = internalInfo[internali].valid(td);

        if (calculate || valid)
        {
            nUnset += !valid;

            distance.primitiveFieldRef()[internali] =
                internalInfo[internali].dist(td);

            (void)std::initializer_list<nil>
            {(
                data.primitiveFieldRef()[internali] =
                    internalInfo[internali].data(td),
                nil()
            ) ... };
        }
    }
    forAll(patchFaceInfo, patchi)
    {
        forAll(patchFaceInfo[patchi], patchFacei)
        {
            const bool valid = patchFaceInfo[patchi][patchFacei].valid(td);

            if (calculate || valid)
            {
                nUnset += !valid;

                distance.boundaryFieldRef()[patchi][patchFacei] =
                    patchFaceInfo[patchi][patchFacei].dist(td) + small;

                (void)std::initializer_list<nil>
                {(
                    data.boundaryFieldRef()[patchi][patchFacei] =
                        patchFaceInfo[patchi][patchFacei].data(td),
                    nil()
                ) ... };
            }
        }
    }

    return nUnset;
}


template<template<class> class PatchField, class GeoMesh>
Foam::label Foam::fvPatchDistWave::calculate
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    GeometricField<scalar, PatchField, GeoMesh>& distance
)
{
    return
        wave<FvWallInfo<wallPoint>>
        (
            mesh,
            getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction),
            -1,
            distance,
            FvFaceCellWave<FvWallInfo<wallPoint>>::defaultTrackingData_
        );
}


template<template<class> class PatchField, class GeoMesh>
void Foam::fvPatchDistWave::correct
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance
)
{
    wave<FvWallInfo<wallFace>>
    (
        mesh,
        getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction),
        nCorrections,
        distance,
        FvFaceCellWave<FvWallInfo<wallFace>>::defaultTrackingData_
    );
}


template<template<class> class PatchField, class GeoMesh>
Foam::label Foam::fvPatchDistWave::calculateAndCorrect
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance
)
{
    const List<labelPair> changedPatchAndFaces =
        getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction);

    const label nUnset =
        wave<FvWallInfo<wallPoint>>
        (
            mesh,
            changedPatchAndFaces,
            -1,
            distance,
            FvFaceCellWave<FvWallInfo<wallPoint>>::defaultTrackingData_
        );

    wave<FvWallInfo<wallFace>>
    (
        mesh,
        changedPatchAndFaces,
        nCorrections,
        distance,
        FvFaceCellWave<FvWallInfo<wallFace>>::defaultTrackingData_
    );

    return nUnset;
}


template
<
    template<class> class WallLocation,
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData
>
Foam::label Foam::fvPatchDistWave::calculate
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td
)
{
    return
        wave<FvWallInfo<WallLocation<wallPoint>>, TrackingData>
        (
            mesh,
            getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction),
            -1,
            distance,
            td,
            data
        );
}


template
<
    template<class> class WallLocation,
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData
>
void Foam::fvPatchDistWave::correct
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td
)
{
    wave<FvWallInfo<WallLocation<wallFace>>, TrackingData>
    (
        mesh,
        getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction),
        nCorrections,
        distance,
        td,
        data
    );
}


template
<
    template<class> class WallLocation,
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData
>
Foam::label Foam::fvPatchDistWave::calculateAndCorrect
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td
)
{
    const List<labelPair> changedPatchAndFaces =
        getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction);

    const label nUnset =
        wave<FvWallInfo<WallLocation<wallPoint>>, TrackingData>
        (
            mesh,
            changedPatchAndFaces,
            -1,
            distance,
            td,
            data
        );

    wave<FvWallInfo<WallLocation<wallFace>>, TrackingData>
    (
        mesh,
        changedPatchAndFaces,
        nCorrections,
        distance,
        td,
        data
    );

    return nUnset;
}


namespace Foam
{
namespace fvPatchDistWave
{
    template<class Type>
    struct WallLocationDataType
    {
        template<class WallLocation>
        using type = WallLocationData<WallLocation, Type>;
    };
}
}


template
<
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData
>
Foam::label Foam::fvPatchDistWave::calculate
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td
)
{
    return
        calculate<WallLocationDataType<DataType>::template type>
        (
            mesh,
            getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction),
            -1,
            distance,
            td,
            data
        );
}


template
<
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData
>
void Foam::fvPatchDistWave::correct
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td
)
{
    correct<WallLocationDataType<DataType>::template type>
    (
        mesh,
        getChangedPatchAndFaces(mesh, patchIDs, minFaceFraction),
        nCorrections,
        distance,
        td,
        data
    );
}


template
<
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData
>
Foam::label Foam::fvPatchDistWave::calculateAndCorrect
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td
)
{
    return
        calculateAndCorrect<WallLocationDataType<DataType>::template type>
        (
            mesh,
            patchIDs,
            minFaceFraction,
            nCorrections,
            distance,
            data,
            td
        );
}


// ************************************************************************* //
