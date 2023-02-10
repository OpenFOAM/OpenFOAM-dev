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

Class
    Foam::fvPatchDistWave

Description
    Takes a set of patches to start FvFaceCellWave from and computed the
    distance at patches and possibly additional transported data.

SourceFiles
    fvPatchDistWave.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatchDistWave_H
#define fvPatchDistWave_H

#include "FvFaceCellWave.H"
#include "volFields.H"
#include "wallPoint.H"
#include "wallFace.H"
#include "WallLocationData.H"
#include "FvWallInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvPatchDistWave
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Get initial set of changed faces
List<labelPair> getChangedPatchAndFaces
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction
);

//- Wave distance (and maybe additional) data from faces. If nCorrections is
//  negative (-1) then the wave propagates through the entire mesh and all
//  values are calculated. If nCorrections is positive, then this many wave
//  steps are computed and the result is corrected only on cells and faces that
//  the wave reaches. Don't use this directly. Use
//  calculate/correct/calculateAndCorrect functions below.
template
<
    class FvWallInfoType,
    class TrackingData,
    template<class> class PatchField,
    class GeoMesh,
    class ... DataType
>
label wave
(
    const fvMesh& mesh,
    const List<labelPair>& changedPatchAndFaces,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    TrackingData& td,
    GeometricField<DataType, PatchField, GeoMesh>& ... data
);

//- Calculate distance data from patches
template<template<class> class PatchField, class GeoMesh>
label calculate
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    GeometricField<scalar, PatchField, GeoMesh>& distance
);

//- Correct distance data from patches
template<template<class> class PatchField, class GeoMesh>
void correct
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance
);

//- Calculate and correct distance data from patches
template<template<class> class PatchField, class GeoMesh>
label calculateAndCorrect
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance
);

//- Calculate distance and additional data from patches, using an
//  arbitrary wall location wave class
template
<
    template<class> class WallLocation,
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData = int
>
label calculate
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td = FvFaceCellWave<nil>::defaultTrackingData_
);

//- Correct distance and additional data from patches, using an
//  arbitrary wall location wave class
template
<
    template<class> class WallLocation,
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData = int
>
void correct
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td = FvFaceCellWave<nil>::defaultTrackingData_
);

//- Calculate and correct distance and additional data from patches, using an
//  arbitrary wall location wave class
template
<
    template<class> class WallLocation,
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData = int
>
label calculateAndCorrect
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td = FvFaceCellWave<nil>::defaultTrackingData_
);

//- Calculate distance and additional data from patches
template
<
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData = int
>
label calculate
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td = FvFaceCellWave<nil>::defaultTrackingData_
);

//- Correct distance and additional data from patches
template
<
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData = int
>
void correct
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td = FvFaceCellWave<nil>::defaultTrackingData_
);

//- Calculate and correct distance and additional data from patches
template
<
    class DataType,
    template<class> class PatchField,
    class GeoMesh,
    class TrackingData = int
>
label calculateAndCorrect
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const scalar minFaceFraction,
    const label nCorrections,
    GeometricField<scalar, PatchField, GeoMesh>& distance,
    GeometricField<DataType, PatchField, GeoMesh>& data,
    TrackingData& td = FvFaceCellWave<nil>::defaultTrackingData_
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvPatchDistWave
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvPatchDistWaveTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
