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

#include "crankConnectingRodMotion.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::Function1s::crankConnectingRodMotion::value
(
    const scalar theta
) const
{
    return
    (
        conRodLength_ + stroke_/2
    )
  - (
        stroke_*cos(degToRad(theta))/2
      + sqrt(sqr(conRodLength_) - sqr(stroke_*sin(degToRad(theta))/2))
    );
}


inline Foam::tmp<Foam::Field<Foam::scalar>>
Foam::Function1s::crankConnectingRodMotion::value
(
    const scalarField&
) const
{
    NotImplemented;
    return tmp<Field<scalar>>(nullptr);
}


inline Foam::scalar Foam::Function1s::crankConnectingRodMotion::integral
(
    const scalar,
    const scalar
) const
{
    NotImplemented;
    return NaN;
}


inline Foam::tmp<Foam::Field<Foam::scalar>>
Foam::Function1s::crankConnectingRodMotion::integral
(
    const scalarField&,
    const scalarField&
) const
{
    NotImplemented;
    return tmp<Field<scalar>>(nullptr);
}


// ************************************************************************* //
