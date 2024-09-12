/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "Repeat.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1s::Repeat<Type>::value(const scalar x) const
{
    const label n = floor((x - start_)/period_);
    const scalar xStar = x - start_ - n*period_;
    return value_->value(xStar);
}


template<class Type>
inline Type Foam::Function1s::Repeat<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    const label n1 = floor((x1 - start_)/period_);
    const label n2 = floor((x2 - start_)/period_);
    const scalar x1Star = x1 - start_ - n1*period_;
    const scalar x2Star = x2 - start_ - n2*period_;
    return
        (n2 - n1)*value_->integral(0, period_)
      + value_->integral(x1Star, x2Star);
}


// ************************************************************************* //
