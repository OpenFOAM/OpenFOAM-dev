/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

#include "Square.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1s::Square<Type>::value(const scalar x) const
{
    const scalar markFraction = markSpace_/(1 + markSpace_);

    const scalar phi = frequency_*(x - start_);
    const scalar fraction = phi - floor(phi);

    return
        amplitude_->value(x)*(fraction < markFraction ? 1 : -1)
      + level_->value(x);
}


template<class Type>
inline Type Foam::Function1s::Square<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    if (!integrable_)
    {
        FatalErrorInFunction
            << "Integration is not defined for " << type() << " functions "
            << "unless the amplitude is constant"
            << exit(FatalError);
    }

    auto integral0 = [&](const scalar x)
    {
        const scalar markFraction = markSpace_/(1 + markSpace_);

        const scalar phi = frequency_*(x - start_);
        const scalar fraction = phi - floor(phi);

        return
            (2*amplitude_->value(x)/frequency_)
           *(
               markFraction*floor(phi + 1 - markFraction)
             + fraction*(fraction < markFraction ? 1 : 0)
             - frequency_/2*x
            );
    };

    return integral0(x2) - integral0(x1) + level_->integral(x1, x2);
}


// ************************************************************************* //
