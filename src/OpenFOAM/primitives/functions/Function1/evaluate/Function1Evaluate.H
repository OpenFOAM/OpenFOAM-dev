/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::Function1Evaluate

Description
    Global functions to evaluate Function1 of GeometricFields

SourceFiles
    Function1Evaluate.C

\*---------------------------------------------------------------------------*/

#ifndef Function1Evaluate_H
#define Function1Evaluate_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void evaluate
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const Function1<Type>& func,
    const GeometricField<Type, PatchField, GeoMesh>& x
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> evaluate
(
    const Function1<Type>& func,
    const dimensionSet& dims,
    const GeometricField<Type, PatchField, GeoMesh>& x
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function1Evaluate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
