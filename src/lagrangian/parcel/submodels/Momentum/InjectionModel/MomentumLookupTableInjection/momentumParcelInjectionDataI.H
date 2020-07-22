/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "momentumParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::point& Foam::momentumParcelInjectionData::x() const
{
    return x_;
}


inline const Foam::vector& Foam::momentumParcelInjectionData::U() const
{
    return U_;
}


inline Foam::scalar Foam::momentumParcelInjectionData::d() const
{
    return d_;
}


inline Foam::scalar Foam::momentumParcelInjectionData::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::momentumParcelInjectionData::mDot() const
{
    return mDot_;
}


inline Foam::point& Foam::momentumParcelInjectionData::x()
{
    return x_;
}


inline Foam::vector& Foam::momentumParcelInjectionData::U()
{
    return U_;
}


inline Foam::scalar& Foam::momentumParcelInjectionData::d()
{
    return d_;
}


inline Foam::scalar& Foam::momentumParcelInjectionData::rho()
{
    return rho_;
}


inline Foam::scalar& Foam::momentumParcelInjectionData::mDot()
{
    return mDot_;
}


// ************************************************************************* //
