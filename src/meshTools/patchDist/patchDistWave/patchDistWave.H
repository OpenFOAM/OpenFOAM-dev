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

Class
    Foam::patchDistWave

Description
    Takes a set of patches to start FaceCellWave from and computed the distance
    at patches and possibly additional transported data.

SourceFiles
    patchDistWave.C

\*---------------------------------------------------------------------------*/

#ifndef patchDistWave_H
#define patchDistWave_H

#include "polyMesh.H"
#include "FaceCellWave.H"
#include "Field.H"
#include "FieldField.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistWave
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Set initial set of changed faces
template<class PatchPointType, class ... InitialPatchData>
void setChangedFaces
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    labelList& changedFaces,
    List<PatchPointType>& faceDist,
    const InitialPatchData& ... initialPatchData
);

//- Copy FaceCellWave values into the cell fields
template
<
    class PatchPointType,
    class TrackingData,
    class DataType,
    class DataMethod
>
label getCellValues
(
    FaceCellWave<PatchPointType, TrackingData>& waveInfo,
    Field<DataType>& cellValues,
    DataMethod method,
    const DataType& stabiliseValue = pTraits<DataType>::zero
);

//- Wave distance data from the patches to the cells
template<class PatchPointType, class TrackingData = int>
label wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance,
    bool correct = true,
    TrackingData& td = FaceCellWave<PatchPointType>::defaultTrackingData_
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistWave
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchDistWaveTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
