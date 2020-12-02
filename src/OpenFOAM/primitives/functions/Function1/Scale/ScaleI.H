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

#include "Scale.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1s::Scale<Type>::value(const scalar x) const
{
    const scalar sx = xScale_->value(x)*x;
    return scale_->value(sx)*value_->value(sx);
}


template<class Type>
inline Type Foam::Function1s::Scale<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    if (integrableScale_)
    {
        const scalar sx = xScale_->value(NaN);
        return scale_->value(NaN)*value_->integral(sx*x1, sx*x2)/sx;
    }

    if (integrableValue_)
    {
        const scalar sx = xScale_->value(NaN);
        return scale_->integral(sx*x1, sx*x2)/sx*value_->value(NaN);
    }

    FatalErrorInFunction
        << "Integration is not defined for " << type() << " functions "
        << "unless x-scaling is constant and either scaling is constant or "
        << "the value is constant"
        << exit(FatalError);

    return Zero;
}


// ************************************************************************* //
