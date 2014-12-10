/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "NoBinaryCollision.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoBinaryCollision<CloudType>::NoBinaryCollision
(
    const dictionary& dict,
    CloudType& cloud
)
:
    BinaryCollisionModel<CloudType>(cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoBinaryCollision<CloudType>::~NoBinaryCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoBinaryCollision<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::NoBinaryCollision<CloudType>::sigmaTcR
(
    const typename CloudType::parcelType& pP,
    const typename CloudType::parcelType& pQ
) const
{
    FatalErrorIn
    (
        "Foam::scalar Foam::NoBinaryCollision<CloudType>::sigmaTcR"
        "("
            "const typename CloudType::parcelType&, "
            "const typename CloudType::parcelType"
        ") const"
    )
        << "sigmaTcR called on NoBinaryCollision model, this should "
        << "not happen, this is not an actual model." << nl
        << "Enclose calls to sigmaTcR within a if (binaryCollision().active()) "
        << " check."
        << abort(FatalError);

    return 0.0;
}


template<class CloudType>
void Foam::NoBinaryCollision<CloudType>::collide
(
    typename CloudType::parcelType& pP,
    typename CloudType::parcelType& pQ
)
{}


// ************************************************************************* //
