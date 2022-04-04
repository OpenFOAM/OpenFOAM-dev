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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvPatchDistWave
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Set initial set of changed faces
template<class PatchPointType, class ... InitialPatchData>
void setChangedFaces
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    List<labelPair>& changedPatchAndFaces,
    List<PatchPointType>& faceDist,
    const InitialPatchData& ... initialPatchData
);

//- Copy FvFaceCellWave values into the cell fields
template<class PatchPointType, class DataType, class DataMethod>
label getCellValues
(
    FvFaceCellWave<PatchPointType>& waveInfo,
    Field<DataType>& cellValues,
    DataMethod method,
    const DataType& stabiliseValue = pTraits<DataType>::zero
);

//- Copy FvFaceCellWave values into the patch field-fields
template<class PatchPointType, class DataType, class DataMethod>
label getPatchValues
(
    FvFaceCellWave<PatchPointType>& wave,
    GeometricBoundaryField<DataType, fvPatchField, volMesh>& valuesBf,
    DataMethod method,
    const DataType& stabiliseValue = pTraits<DataType>::zero
);

//- Wave distance data from the patches to the cells and other patch faces
template<class PatchPointType>
label wave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    volScalarField& distance,
    bool correct = true
);

//- Wave distance and auxiliary data from the patches to the cells and other
//  patch faces
template<class PatchPointType>
label wave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const GeometricBoundaryField
        <typename PatchPointType::dataType, fvPatchField, volMesh>&
        initialPatchData,
    volScalarField& distance,
    VolField<typename PatchPointType::dataType>& data,
    bool correct = true
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
