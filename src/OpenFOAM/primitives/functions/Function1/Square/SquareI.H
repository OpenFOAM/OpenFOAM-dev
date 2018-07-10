/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1Types::Square<Type>::value(const scalar t) const
{
    // Number of waves including fractions
    scalar waves = frequency_->value(t)*(t - t0_);

    // Number of complete waves
    scalar nWaves;

    // Fraction of last incomplete wave
    scalar waveFrac = std::modf(waves, &nWaves);

    // Mark fraction of a wave
    scalar markFrac = markSpace_/(1.0 + markSpace_);

    return
        amplitude_->value(t)
       *(waveFrac < markFrac ? 1 : -1)
       *scale_->value(t)
      + level_->value(t);
}


// ************************************************************************* //
