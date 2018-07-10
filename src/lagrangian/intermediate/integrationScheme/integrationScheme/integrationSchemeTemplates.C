/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "integrationScheme.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::integrationScheme::explicitDelta
(
    const Type& phi,
    const scalar dtEff,
    const Type& Alpha,
    const scalar Beta
)
{
    return (Alpha - Beta*phi)*dtEff;
}


template<class Type>
inline Type Foam::integrationScheme::delta
(
    const Type& phi,
    const scalar dt,
    const Type& Alpha,
    const scalar Beta
) const
{
    return explicitDelta(phi, dtEff(dt, Beta), Alpha, Beta);
}


template<class Type>
inline Type Foam::integrationScheme::partialDelta
(
    const Type& phi,
    const scalar dt,
    const Type& Alpha,
    const scalar Beta,
    const Type& alphai,
    const scalar betai
) const
{
    return
        explicitDelta(phi, dt, alphai, betai)
      - explicitDelta(phi, sumDtEff(dt, Beta), Alpha, Beta)*betai;
}


// ************************************************************************* //
