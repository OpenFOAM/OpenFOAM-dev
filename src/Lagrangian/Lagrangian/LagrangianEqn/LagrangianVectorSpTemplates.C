/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianSp.H"
#include "LagrangianEqn.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<template<class> class PrimitiveField>
Foam::LagrangianSp<Foam::vector>::LagrangianSp
(
    const LagrangianEqnBase& eqn,
    const LagrangianSubField<scalar, PrimitiveField>& S
)
:
    scalarCoeff_(eqn, S),
    tensorCoeff_(eqn)
{}


template<template<class> class PrimitiveField>
Foam::LagrangianSp<Foam::vector>::LagrangianSp
(
    const LagrangianEqnBase& eqn,
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
:
    scalarCoeff_(eqn, tS),
    tensorCoeff_(eqn)
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator+=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    if (tensorCoeff_.valid())
    {
        tensorCoeff_ += tensor::I*S;
    }
    else
    {
        scalarCoeff_ += S;
    }
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator+=
(
    const LagrangianSubField<tensor, PrimitiveField>& S
)
{
    makeTensor();
    tensorCoeff_ += S;
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator+=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    if (tensorCoeff_.valid())
    {
        tensorCoeff_ += tensor::I*tS;
    }
    else
    {
        scalarCoeff_ += tS;
    }
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator+=
(
    const tmp<LagrangianSubField<tensor, PrimitiveField>>& tS
)
{
    makeTensor();
    tensorCoeff_ += tS;
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator-=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    if (tensorCoeff_.valid())
    {
        tensorCoeff_ -= tensor::I*S;
    }
    else
    {
        scalarCoeff_ -= S;
    }
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator-=
(
    const LagrangianSubField<tensor, PrimitiveField>& S
)
{
    makeTensor();
    tensorCoeff_ -= S;
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator-=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    if (tensorCoeff_.valid())
    {
        tensorCoeff_ -= tensor::I*tS;
    }
    else
    {
        scalarCoeff_ -= tS;
    }
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator-=
(
    const tmp<LagrangianSubField<tensor, PrimitiveField>>& tS
)
{
    makeTensor();
    tensorCoeff_ -= tS;
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator*=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    scalarCoeff_ *= S;
    tensorCoeff_ *= S;
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator*=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    *this *= tS();
    tS.clear();
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator/=
(
    const LagrangianSubField<scalar, PrimitiveField>& S
)
{
    scalarCoeff_ /= S;
    tensorCoeff_ /= S;
}


template<template<class> class PrimitiveField>
void Foam::LagrangianSp<Foam::vector>::operator/=
(
    const tmp<LagrangianSubField<scalar, PrimitiveField>>& tS
)
{
    *this /= tS();
    tS.clear();
}


// ************************************************************************* //
