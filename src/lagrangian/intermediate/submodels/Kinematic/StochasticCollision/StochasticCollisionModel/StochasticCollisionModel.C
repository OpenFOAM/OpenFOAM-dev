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

#include "StochasticCollisionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::StochasticCollisionModel<CloudType>::collide(const scalar dt)
{
    notImplemented
    (
        "void Foam::NoStochasticCollision<CloudType>::collide(const scalar)"
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::StochasticCollisionModel
(
    CloudType& owner
)
:
    CloudSubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::StochasticCollisionModel
(
    const StochasticCollisionModel<CloudType>& cm
)
:
    CloudSubModelBase<CloudType>(cm)
{}


template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::StochasticCollisionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticCollisionModel<CloudType>::~StochasticCollisionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::StochasticCollisionModel<CloudType>::update(const scalar dt)
{
    if (this->active())
    {
        this->collide(dt);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "StochasticCollisionModelNew.C"

// ************************************************************************* //
