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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::fv::cellSetOption::timeStart() const
{
    return timeStart_;
}


inline Foam::scalar Foam::fv::cellSetOption::duration() const
{
    return duration_;
}


inline bool Foam::fv::cellSetOption::inTimeLimits(const scalar time) const
{
    return
    (
        (timeStart_ < 0)
     ||
        (
            (mesh_.time().value() >= timeStart_)
         && (mesh_.time().value() <= (timeStart_ + duration_))
        )
    );
}


inline const Foam::fv::cellSetOption::selectionModeType&
Foam::fv::cellSetOption::selectionMode() const
{
    return selectionMode_;
}


inline const Foam::word& Foam::fv::cellSetOption::cellSetName() const
{
    return cellSetName_;
}


inline Foam::scalar Foam::fv::cellSetOption::V() const
{
    return V_;
}


inline const Foam::labelList& Foam::fv::cellSetOption::cells() const
{
    return cells_;
}


inline Foam::scalar& Foam::fv::cellSetOption::timeStart()
{
    return timeStart_;
}


inline Foam::scalar& Foam::fv::cellSetOption::duration()
{
    return duration_;
}


// ************************************************************************* //
