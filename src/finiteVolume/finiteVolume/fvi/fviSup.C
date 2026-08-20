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

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvi
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolInternalField<Type>>
Su
(
    const VolInternalField<Type>& su,
    const VolField<Type>& vf
)
{
    return su;
}

template<class Type, class Expression, class>
tmp<VolInternalField<Type>> Su
(
    const Expression& su,
    const VolField<Type>&
)
{
    return eval(su);
}

template<class Type>
tmp<VolInternalField<Type>>
Su
(
    const tmp<VolInternalField<Type>>& tsu,
    const VolField<Type>& vf
)
{
    return tsu;
}


template<class Type>
tmp<VolInternalField<Type>>
Sp
(
    const VolInternalField<scalar>& sp,
    const VolField<Type>& vf
)
{
    return sp*vf;
}

template<class Type, class Expression, class>
tmp<VolInternalField<Type>> Sp
(
    const Expression& sp,
    const VolField<Type>& vf
)
{
    return eval(sp*vf.internalField());
}

template<class Type>
tmp<VolInternalField<Type>>
Sp
(
    const tmp<VolInternalField<scalar>>& tsp,
    const VolField<Type>& vf
)
{
    return tsp*vf();
}


template<class Type>
tmp<VolInternalField<Type>>
Sp
(
    const dimensionedScalar& sp,
    const VolField<Type>& vf
)
{
    return sp*vf;
}


template<class Type>
tmp<VolInternalField<Type>>
SuSp
(
    const VolInternalField<scalar>& sp,
    const VolField<Type>& vf
)
{
    return sp*vf;
}

template<class Type, class Expression, class>
tmp<VolInternalField<Type>> SuSp
(
    const Expression& sp,
    const VolField<Type>& vf
)
{
    return eval(sp*vf.internalField());
}

template<class Type>
tmp<VolInternalField<Type>>
SuSp
(
    const tmp<VolInternalField<scalar>>& tsp,
    const VolField<Type>& vf
)
{
    return tsp*vf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
