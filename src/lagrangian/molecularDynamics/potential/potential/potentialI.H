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

inline Foam::label Foam::potential::nIds() const
{
    return idList_.size();
}


inline const Foam::List<Foam::word>& Foam::potential::idList() const
{
    return idList_;
}


inline const Foam::List<Foam::word>& Foam::potential::siteIdList() const
{
    return siteIdList_;
}


inline Foam::scalar Foam::potential::potentialEnergyLimit() const
{
    return potentialEnergyLimit_;
}


inline Foam::label Foam::potential::nPairPotentials() const
{
    return pairPotentials_.size();
}


inline const Foam::labelList& Foam::potential::removalOrder() const
{
    return removalOrder_;
}


inline const Foam::pairPotentialList& Foam::potential::pairPotentials() const
{
    return pairPotentials_;
}


inline const Foam::tetherPotentialList&
Foam::potential::tetherPotentials() const
{
    return tetherPotentials_;
}


inline const Foam::vector& Foam::potential::gravity() const
{
    return gravity_;
}


// ************************************************************************* //
