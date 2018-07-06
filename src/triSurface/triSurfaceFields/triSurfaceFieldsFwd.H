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

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceFieldsFwd_H
#define triSurfaceFieldsFwd_H

#include "fieldTypes.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
class DimensionedField;

class triSurfaceGeoMesh;

typedef Foam::DimensionedField<label, triSurfaceGeoMesh>
    triSurfaceLabelField;
typedef Foam::DimensionedField<scalar, triSurfaceGeoMesh>
    triSurfaceScalarField;
typedef Foam::DimensionedField<vector, triSurfaceGeoMesh>
    triSurfaceVectorField;
typedef Foam::DimensionedField<sphericalTensor, triSurfaceGeoMesh>
    triSurfaceSphericalTensorField;
typedef Foam::DimensionedField<symmTensor, triSurfaceGeoMesh>
    triSurfaceSymmTensorField;
typedef Foam::DimensionedField<tensor, triSurfaceGeoMesh>
    triSurfaceTensorField;

class triSurfacePointGeoMesh;

typedef Foam::DimensionedField<label, triSurfacePointGeoMesh>
    triSurfacePointLabelField;
typedef Foam::DimensionedField<scalar, triSurfacePointGeoMesh>
    triSurfacePointScalarField;
typedef Foam::DimensionedField<vector, triSurfacePointGeoMesh>
    triSurfacePointVectorField;
typedef Foam::DimensionedField<sphericalTensor, triSurfacePointGeoMesh>
    triSurfacePointSphericalTensorField;
typedef Foam::DimensionedField<symmTensor, triSurfacePointGeoMesh>
    triSurfacePointSymmTensorField;
typedef Foam::DimensionedField<tensor, triSurfacePointGeoMesh>
    triSurfacePointTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
