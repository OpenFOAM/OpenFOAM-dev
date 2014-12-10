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

#include "NoPhaseChange.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoPhaseChange<CloudType>::NoPhaseChange
(
    const dictionary&,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoPhaseChange<CloudType>::NoPhaseChange
(
    const NoPhaseChange<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoPhaseChange<CloudType>::~NoPhaseChange()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoPhaseChange<CloudType>::active() const
{
    return false;
}


template<class CloudType>
void Foam::NoPhaseChange<CloudType>::calculate
(
    const scalar,
    const label,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    scalarField&
) const
{
    // Nothing to do...
}


// ************************************************************************* //
