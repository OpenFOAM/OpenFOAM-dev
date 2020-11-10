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

\*---------------------------------------------------------------------------*/

#include "Function1Evaluate.H"

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::evaluate
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const Function1<Type>& func,
    const GeometricField<Type, PatchField, GeoMesh>& f
)
{
    result.primitiveFieldRef() = func.value(f());

    typename GeometricField<Type, PatchField, GeoMesh>::Boundary& bresult =
        result.boundaryFieldRef();

    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary& bf =
        f.boundaryField();

    forAll(bresult, patchi)
    {
        bresult[patchi] = func.value(bf[patchi]);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>> Foam::evaluate
(
    const Function1<Type>& func,
    const dimensionSet& dims,
    const GeometricField<Type, PatchField, GeoMesh>& f
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tresult
    (
        GeometricField<Type, PatchField, GeoMesh>::New
        (
            func.name() + "(" + f.name() + ')',
            f.mesh(),
            dims
        )
    );

    evaluate(tresult.ref(), func, f);

    return tresult;
}


// ************************************************************************* //
