/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "DimensionedTensorField.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<tensor, GeoMesh, Field>> inv
(
    const DimensionedField<tensor, GeoMesh, PrimitiveField>& dtf,
    const Vector<label>& solutionD
)
{
    tmp<DimensionedField<tensor, GeoMesh, Field>> tRes
    (
        DimensionedField<tensor, GeoMesh, Field>::New
        (
            "inv(" + dtf.name() + ',' + name(solutionD) + ')',
            dtf.mesh(),
            inv(dtf.dimensions())
        )
    );

    inv(tRes.ref().primitiveFieldRef(), dtf.primitiveField(), solutionD);

    return tRes;
}


template<class GeoMesh, template<class> class PrimitiveField>
tmp<DimensionedField<tensor, GeoMesh, Field>> inv
(
    const tmp<DimensionedField<tensor, GeoMesh, PrimitiveField>>& tdtf,
    const Vector<label>& solutionD
)
{
    tmp<DimensionedField<tensor, GeoMesh, Field>> tRes(inv(tdtf(), solutionD));
    tdtf.clear();
    return tRes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
