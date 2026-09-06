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

InClass
    Foam::DimensionedSymmTensorField

Description
    SymmTensor specific part of the implementation of DimensionedField.

SourceFiles
    DimensionedSymmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedSymmTensorField_H
#define DimensionedSymmTensorField_H

#include "DimensionedTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<symmTensor, GeoMesh, Field>> inv
(
    const DimensionedField<symmTensor, GeoMesh, PrimitiveField>&,
    const Vector<label>& solutionD
);

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<symmTensor, GeoMesh, Field>> inv
(
    const tmp<DimensionedField<symmTensor, GeoMesh, PrimitiveField>>& tdtf,
    const Vector<label>& solutionD
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DimensionedSymmTensorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
