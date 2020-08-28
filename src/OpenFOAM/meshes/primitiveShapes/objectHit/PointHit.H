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

Class
    Foam::PointHit

Description
    This class describes the interaction of a face and a point. It
    carries the info of a successful hit and (if successful), returns
    the interaction point.

\*---------------------------------------------------------------------------*/

#ifndef PointHit_H
#define PointHit_H

#include "bool.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class Ostream;


// Forward declaration of friend functions and operators

template<class Point> class PointHit;

template<class Point>
inline Ostream& operator<<(Ostream&, const PointHit<Point>&);


/*---------------------------------------------------------------------------*\
                          Class PointHit Declaration
\*---------------------------------------------------------------------------*/

template<class Point>
class PointHit
{
    // Private Data

        //- Hit success
        bool hit_;

        //- Point of hit; for miss holds best estimate outside the object
        Point hitPoint_;

        //- Distance to hit point
        scalar distance_;

        //- Eligible miss
        bool eligibleMiss_;


public:

    // Constructors

        //- Construct null
        PointHit()
        :
            hit_(false),
            hitPoint_(Zero),
            distance_(great),
            eligibleMiss_(false)
        {}

        //- Construct from components
        PointHit
        (
            const bool hit,
            const Point& p,
            const scalar dist,
            const bool eligibleMiss
        )
        :
            hit_(hit),
            hitPoint_(p),
            distance_(dist),
            eligibleMiss_(eligibleMiss)
        {}

        //- Construct from point. Hit and distance set later
        PointHit(const Point& p)
        :
            hit_(false),
            hitPoint_(p),
            distance_(great),
            eligibleMiss_(false)
        {}


    // Member Functions

        //- Is there a hit
        bool hit() const
        {
            return hit_;
        }

        //- Return hit point
        const Point& hitPoint() const
        {
            if (!hit_)
            {
                FatalErrorInFunction
                    << "requested a hit point for a miss"
                    << abort(FatalError);
            }

            return hitPoint_;
        }

        //- Return distance to hit
        scalar distance() const
        {
            return distance_;
        }

        //- Return miss point
        const Point& missPoint() const
        {
            if (hit_)
            {
                FatalErrorInFunction
                    << "requested a miss point for a hit"
                    << abort(FatalError);
            }

            return hitPoint_;
        }

        //- Return point with no checking
        const Point& rawPoint() const
        {
            return hitPoint_;
        }

        //- Is this an eligible miss
        bool eligibleMiss() const
        {
            return eligibleMiss_;
        }

        void setHit()
        {
            hit_ = true;
            eligibleMiss_ = false;
        }

        void setMiss(const bool eligible)
        {
            hit_ = false;
            eligibleMiss_ = eligible;
        }

        void setPoint(const Point& p)
        {
            hitPoint_ = p;
        }

        void setDistance(const scalar d)
        {
            distance_ = d;
        }


    // Ostream operator

        friend Ostream& operator<< <Point>
        (
            Ostream& os,
            const PointHit<Point>& b
        );
};


template<class Point>
inline Ostream& operator<<(Ostream& os, const PointHit<Point>& b)
{
    os  << b.hit() << token::SPACE
        << b.rawPoint() << token::SPACE
        << b.distance() << token::SPACE
        << b.eligibleMiss();

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
