/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "cellPoint_interpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolations::cellPointBase<Type>::cellPointBase
(
    const VolField<Type>& psi
)
:
    fieldInterpolation<Type, cellPoint<Type>>(psi),
    volPointInterpolation<Type>(psi)
{
    // The interpolation process needs the tetrahedral decomposition
    (void)psi.mesh().tetBasePtIs();
}


template<class Type>
Foam::interpolations::cellPointBase<Type>::cellPointBase
(
    const cellPointBase<Type>& i
)
:
    fieldInterpolation<Type, cellPoint<Type>>(i),
    volPointInterpolation<Type>(i)
{
    // The interpolation process needs the tetrahedral decomposition
    (void)this->psi().mesh().tetBasePtIs();
}


template<class Type>
Foam::interpolations::cellPointBase<Type>::cellPointBase
(
    const VolField<Type>& psi,
    tmp<PointField<Type>> psip
)
:
    fieldInterpolation<Type, cellPoint<Type>>(psi),
    volPointInterpolation<Type>(psi, psip)
{
    // The interpolation process needs the tetrahedral decomposition
    (void)psi.mesh().tetBasePtIs();
}


// ************************************************************************* //
