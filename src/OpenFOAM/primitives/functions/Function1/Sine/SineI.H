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

#include "Sine.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1s::Sine<Type>::value(const scalar x) const
{
    const scalar phi = constant::mathematical::twoPi*frequency_*(x - start_);

    return amplitude_->value(x)*sin(phi) + level_->value(x);
}


template<class Type>
inline Type Foam::Function1s::Sine<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    using namespace constant::mathematical;

    if (!integrable_)
    {
        FatalErrorInFunction
            << "Integration is not defined for " << type() << " functions "
            << "unless the amplitude is constant"
            << exit(FatalError);
    }

    const scalar phi1 = constant::mathematical::twoPi*frequency_*(x1 - start_);
    const scalar phi2 = constant::mathematical::twoPi*frequency_*(x2 - start_);

    return
      - amplitude_->value(NaN)*(cos(phi2) - cos(phi1))/(twoPi*frequency_)
      + level_->integral(x1, x2);
}


// ************************************************************************* //
