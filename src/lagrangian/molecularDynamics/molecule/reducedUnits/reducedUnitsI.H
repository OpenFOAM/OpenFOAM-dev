/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

inline Foam::scalar Foam::reducedUnits::refLength() const
{
    return refLength_;
}


inline Foam::scalar Foam::reducedUnits::refTime() const
{
    return refTime_;
}


inline Foam::scalar Foam::reducedUnits::refMass() const
{
    return refMass_;
}


inline Foam::scalar Foam::reducedUnits::refTemp() const
{
    return refTemp_;
}


inline Foam::scalar Foam::reducedUnits::refEnergy() const
{
    return refEnergy_;
}


inline Foam::scalar Foam::reducedUnits::refForce() const
{
    return refForce_;
}


inline Foam::scalar Foam::reducedUnits::refVelocity() const
{
    return refVelocity_;
}


inline Foam::scalar Foam::reducedUnits::refVolume() const
{
    return refVolume_;
}


inline Foam::scalar Foam::reducedUnits::refPressure() const
{
    return refPressure_;
}


inline Foam::scalar Foam::reducedUnits::refMassDensity() const
{
    return refMassDensity_;
}


inline Foam::scalar Foam::reducedUnits::refNumberDensity() const
{
    return refNumberDensity_;
}


// ************************************************************************* //
