/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "Spline.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Derived>
inline const Foam::pointField& Foam::Spline<Derived>::knots() const
{
    return knots_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Derived>
inline Foam::label Foam::Spline<Derived>::nSegments() const
{
    return knots_.size() - 1;
}


template<class Derived>
inline const Foam::point& Foam::Spline<Derived>::start() const
{
    return knots_.first();
}


template<class Derived>
inline const Foam::point& Foam::Spline<Derived>::end() const
{
    return knots_.last();
}


// ************************************************************************* //
