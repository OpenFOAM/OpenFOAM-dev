/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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

#include "AntoineExtended.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::Function1s::AntoineExtended::value
(
    const scalar T
) const
{
    return exp(A_ + B_/(C_ + T) + D_*T + E_*log(T) + F_*pow(T, G_));
}


inline Foam::scalar Foam::Function1s::AntoineExtended::derivative
(
    const scalar T
) const
{
    return value(T)*(- B_/sqr(C_ + T) + D_ + E_/T + F_*G_*pow(T, G_ - 1));
}


inline Foam::scalar Foam::Function1s::AntoineExtended::integral
(
    const scalar T1,
    const scalar T2
) const
{
    NotImplemented;
    return NaN;
}


// ************************************************************************* //
