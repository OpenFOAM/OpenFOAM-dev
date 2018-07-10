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

#include "TrajectoryCollision.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::TrajectoryCollision<CloudType>::collide
(
    typename CloudType::parcelType::trackingData& td,
    const scalar dt
)
{
    ORourkeCollision<CloudType>::collide(td, dt);
}


template<class CloudType>
bool Foam::TrajectoryCollision<CloudType>::collideParcels
(
    const scalar dt,
    parcelType& p1,
    parcelType& p2,
    scalar& m1,
    scalar& m2
)
{
    bool coalescence = false;

    const vector& pos1 = p1.position();
    const vector& pos2 = p2.position();

    const vector& U1 = p1.U();
    const vector& U2 = p2.U();

    vector URel = U1 - U2;

    vector d = pos2 - pos1;
    scalar magd = mag(d);

    scalar vAlign = URel & (d/(magd + rootVSmall));

    if (vAlign > 0)
    {
        const scalar d1 = p1.d();
        const scalar d2 = p2.d();

        scalar sumD = d1 + d2;

        if (vAlign*dt > magd - 0.5*sumD)
        {
            scalar magU1 = mag(U1) + rootVSmall;
            scalar magU2 = mag(U2) + rootVSmall;
            vector n1 = U1/magU1;
            vector n2 = U2/magU2;

            scalar n1n2 = n1 & n2;
            scalar n1d = n1 & d;
            scalar n2d = n2 & d;

            scalar det = 1.0 - sqr(n1n2);

            scalar alpha = great;
            scalar beta = great;

            if (mag(det) > 1.0e-4)
            {
                beta = -(n2d - n1n2*n1d)/det;
                alpha = n1d + n1n2*beta;
            }

            alpha /= magU1*dt;
            beta /= magU2*dt;

            // is collision possible within this timestep
            if ((alpha > 0) && (alpha < 1.0) && (beta > 0) && (beta < 1.0))
            {
                vector p1c = pos1 + alpha*U1*dt;
                vector p2c = pos2 + beta*U2*dt;

                scalar closestDist = mag(p1c - p2c);

                scalar collProb =
                    pow(0.5*sumD/max(0.5*sumD, closestDist), cSpace_)
                   *exp(-cTime_*mag(alpha - beta));

                scalar xx = this->owner().rndGen().template sample01<scalar>();

                // collision occurs
                if (xx > collProb)
                {
                    if (d1 > d2)
                    {
                        coalescence = this->collideSorted(dt, p1, p2, m1, m2);
                    }
                    else
                    {
                        coalescence = this->collideSorted(dt, p2, p1, m2, m1);
                    }
                }
            }
        }
    }

    return coalescence;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TrajectoryCollision<CloudType>::TrajectoryCollision
(
    const dictionary& dict,
    CloudType& owner
)
:
    ORourkeCollision<CloudType>(dict, owner, typeName),
    cSpace_(readScalar(this->coeffDict().lookup("cSpace"))),
    cTime_(readScalar(this->coeffDict().lookup("cTime")))
{}


template<class CloudType>
Foam::TrajectoryCollision<CloudType>::TrajectoryCollision
(
    const TrajectoryCollision<CloudType>& cm
)
:
    ORourkeCollision<CloudType>(cm),
    cSpace_(cm.cSpace_),
    cTime_(cm.cTime_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TrajectoryCollision<CloudType>::~TrajectoryCollision()
{}


// ************************************************************************* //
