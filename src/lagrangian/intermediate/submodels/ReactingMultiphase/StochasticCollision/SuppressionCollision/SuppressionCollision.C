/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "SuppressionCollision.H"
#include "kinematicCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SuppressionCollision<CloudType>::collide(const scalar dt)
{
    const kinematicCloud& sc =
        this->owner().mesh().template
        lookupObject<kinematicCloud>(suppressionCloud_);

    volScalarField vDotSweep(sc.vDotSweep());

    dimensionedScalar Dt("dt", dimTime, dt);
    volScalarField P(type() + ":p", 1.0 - exp(-vDotSweep*Dt));

    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();
        label cellI = p.cell();

        scalar xx = this->owner().rndGen().template sample01<scalar>();

        if (xx < P[cellI])
        {
            p.canCombust() = -1;
            p.typeId() = max(p.typeId(), suppressedParcelType_);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SuppressionCollision<CloudType>::SuppressionCollision
(
    const dictionary& dict,
    CloudType& owner
)
:
    StochasticCollisionModel<CloudType>(dict, owner, typeName),
    suppressionCloud_(this->coeffDict().lookup("suppressionCloud")),
    suppressedParcelType_
    (
        this->coeffDict().lookupOrDefault("suppressedParcelType", -1)
    )
{}


template<class CloudType>
Foam::SuppressionCollision<CloudType>::SuppressionCollision
(
    const SuppressionCollision<CloudType>& cm
)
:
    StochasticCollisionModel<CloudType>(cm),
    suppressionCloud_(cm.suppressionCloud_),
    suppressedParcelType_(cm.suppressedParcelType_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SuppressionCollision<CloudType>::~SuppressionCollision()
{}


// ************************************************************************* //
