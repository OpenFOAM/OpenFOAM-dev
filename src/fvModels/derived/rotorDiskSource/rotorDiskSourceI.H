/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "rotorDiskSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::fv::rotorDiskSource::rhoRef() const
{
    return rhoRef_;
}


inline Foam::scalar Foam::fv::rotorDiskSource::omega() const
{
    return omega_;
}


inline const Foam::fvCellSet& Foam::fv::rotorDiskSource::set() const
{
    return set_;
}


inline const Foam::List<Foam::point>& Foam::fv::rotorDiskSource::x() const
{
    return x_;
}


inline const Foam::coordinateSystems::cylindrical&
Foam::fv::rotorDiskSource::coordSys() const
{
    return coordSys_;
}


// ************************************************************************* //
