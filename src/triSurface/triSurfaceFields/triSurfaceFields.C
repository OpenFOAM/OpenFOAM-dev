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

#include "triSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const word triSurfaceLabelField::typeName("triSurfaceLabelField");

template<>
const word triSurfaceScalarField::typeName("triSurfaceScalarField");

template<>
const word triSurfaceVectorField::typeName("triSurfaceVectorField");

template<>
const word triSurfaceSphericalTensorField::typeName
("triSurfaceSphericalTensorField");

template<>
const word triSurfaceSymmTensorField::typeName
("triSurfaceSymmTensorField");

template<>
const word triSurfaceTensorField::typeName("triSurfaceTensorField");


template<>
const word triSurfacePointLabelField::typeName("triSurfacePointLabelField");

template<>
const word triSurfacePointScalarField::typeName("triSurfacePointScalarField");

template<>
const word triSurfacePointVectorField::typeName("triSurfacePointVectorField");

template<>
const word triSurfacePointSphericalTensorField::typeName
("triSurfacePointSphericalTensorField");

template<>
const word triSurfacePointSymmTensorField::typeName
("triSurfacePointSymmTensorField");

template<>
const word triSurfacePointTensorField::typeName("triSurfacePointTensorField");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
