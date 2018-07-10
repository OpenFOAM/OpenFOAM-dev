/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::slicedVolFields

\*---------------------------------------------------------------------------*/

#ifndef slicedVolFieldsFwd_H
#define slicedVolFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class volMesh;

template<class Type>
class fvPatchField;

template<class Type>
class slicedFvPatchField;

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
class SlicedGeometricField;

typedef
SlicedGeometricField<scalar, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolScalarField;

typedef
SlicedGeometricField<vector, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolVectorField;

typedef
SlicedGeometricField<sphericalTensor, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolSphericalTensorField;

typedef
SlicedGeometricField<symmTensor, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolSymmTensorField;

typedef
SlicedGeometricField<tensor, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
