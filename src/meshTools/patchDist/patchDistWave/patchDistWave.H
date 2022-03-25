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

#include "Field.H"
#include "FieldField.H"
#include "UPtrList.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class wallPoint;

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
template<class PatchPointType, class DataType, class DataMethod>
label getCellValues
(
    const polyMesh& mesh,
    FaceCellWave<PatchPointType>& waveInfo,
    Field<DataType>& cellValues,
    DataMethod method,
    const DataType& stabiliseValue = pTraits<DataType>::zero
);

//- Copy FaceCellWave values into the patch field-fields
template
<
    class PatchPointType,
    template<class> class PatchField,
    class DataType,
    class DataMethod
>
label getPatchValues
(
    const polyMesh& mesh,
    FaceCellWave<PatchPointType>& wave,
    FieldField<PatchField, DataType>& patchValues,
    DataMethod method,
    const DataType& stabiliseValue = pTraits<DataType>::zero
);

//- Wave distance data from the patches to the cells
template<class PatchPointType>
label wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance,
    bool correct = true
);

//- Wave distance and auxiliary data from the patches to the cells and other
//  patch faces
template<class PatchPointType, template<class> class PatchField>
label wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    const FieldField<PatchField, typename PatchPointType::dataType>&
        initialPatchData,
    scalarField& cellDistance,
    FieldField<PatchField, scalar>& patchDistance,
    Field<typename PatchPointType::dataType>& cellData,
    FieldField<PatchField, typename PatchPointType::dataType>& patchData,
    bool correct = true
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
