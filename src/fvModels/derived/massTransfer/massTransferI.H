/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "massTransfer.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline Foam::label Foam::fv::massTransfer::index
(
    const Pair<word>& names,
    const word& name
)
{
    return name == names.first() ? 0 : name == names.second() ? 1 : -1;
}


inline Foam::label Foam::fv::massTransfer::sign
(
    const Pair<word>& names,
    const word& name
)
{
    return name == names.first() ? -1 : name == names.second() ? +1 : 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::Pair<Foam::word>& Foam::fv::massTransfer::phaseNames() const
{
    return phaseNames_;
}


inline const Foam::Pair<Foam::word>& Foam::fv::massTransfer::alphaNames() const
{
    return alphaNames_;
}


inline const Foam::Pair<Foam::word>& Foam::fv::massTransfer::rhoNames() const
{
    return rhoNames_;
}


// ************************************************************************* //

