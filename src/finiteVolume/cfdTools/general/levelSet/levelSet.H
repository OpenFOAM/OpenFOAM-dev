/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

Description:
    Functions to create an averaged field from a discontinuous field defined by
    a level-set.

SourceFiles:
    levelSet.C
    levelSetTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef levelSet_H
#define levelSet_H

#include "DimensionedField.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Calculate the average value of two fields, one on each side of a level set
//  The level set and the fields are given both on the points and cell-centres.
template<class Type>
tmp<DimensionedField<Type, volMesh>> levelSetAverage
(
    const fvMesh& mesh,
    const scalarField& levelC,
    const scalarField& levelP,
    const DimensionedField<Type, volMesh>& positiveC,
    const DimensionedField<Type, pointMesh>& positiveP,
    const DimensionedField<Type, volMesh>& negativeC,
    const DimensionedField<Type, pointMesh>& negativeP
);

//- As the above overload, but on the faces of a patch
template<class Type>
tmp<Field<Type>> levelSetAverage
(
    const fvPatch& patch,
    const scalarField& levelF,
    const scalarField& levelP,
    const Field<Type>& positiveF,
    const Field<Type>& positiveP,
    const Field<Type>& negativeF,
    const Field<Type>& negativeP
);

//- Calculate the volume-fraction that a level set occupies. This gives the the
//  same result as levelSetAverage if the fields passed to the latter are
//  uniformly 0 and 1. The above flag flips the direction.
tmp<DimensionedField<scalar, volMesh>> levelSetFraction
(
    const fvMesh& mesh,
    const scalarField& levelC,
    const scalarField& levelP,
    const bool above
);

//- As the above oveload, but on the faces of a patch
tmp<scalarField> levelSetFraction
(
    const fvPatch& patch,
    const scalarField& levelF,
    const scalarField& levelP,
    const bool above
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "levelSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
