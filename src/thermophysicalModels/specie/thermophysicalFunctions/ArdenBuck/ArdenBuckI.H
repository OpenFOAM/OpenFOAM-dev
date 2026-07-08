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

#include "ArdenBuck.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar
Foam::Function1s::ArdenBuck::xByTC(const scalar TC) const
{
    return (B_ - TC/C_)/(D_ + TC);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::Function1s::ArdenBuck::value(const scalar T) const
{
    const scalar TC = T - zeroC_;

    return A_*exp(TC*xByTC(TC));
}


inline Foam::scalar Foam::Function1s::ArdenBuck::derivative
(
    const scalar T
) const
{
    const scalar TC = T - zeroC_;
    const scalar x = xByTC(TC);

    return A_*exp(TC*x)*(D_*x - TC/C_)/(D_ + TC);
}


inline Foam::scalar Foam::Function1s::ArdenBuck::integral
(
    const scalar T1,
    const scalar T2
) const
{
    NotImplemented;
    return NaN;
}


// ************************************************************************* //
