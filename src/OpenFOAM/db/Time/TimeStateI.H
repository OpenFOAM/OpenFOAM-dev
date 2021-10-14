/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::TimeState::timeIndex() const
{
    return timeIndex_;
}


inline Foam::scalar Foam::TimeState::deltaTValue() const
{
    return deltaT_;
}


inline Foam::scalar Foam::TimeState::deltaT0Value() const
{
    return deltaT0_;
}


inline Foam::dimensionedScalar Foam::TimeState::deltaT() const
{
    return dimensionedScalar("deltaT", dimTime, deltaT_);
}


inline Foam::dimensionedScalar Foam::TimeState::deltaT0() const
{
    return dimensionedScalar("deltaT0", dimTime, deltaT0_);
}


inline bool Foam::TimeState::writeTime() const
{
    return writeTime_;
}


inline bool Foam::TimeState::outputTime() const
{
    return writeTime_;
}


// ************************************************************************* //
