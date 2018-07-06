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

Class
    Foam::PointIndexHit

Description
    This class describes the interaction of (usually) a face and a point.
    It carries the info of a successful hit and (if successful),
    returns the interaction point.

    like pointHit but carries face (or cell, edge etc.) index

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef PointIndexHit_H
#define PointIndexHit_H

#include "bool.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PointIndexHit Declaration
\*---------------------------------------------------------------------------*/

template<class Point>
class PointIndexHit
{
    // Private data

        //- Hit success
        bool hit_;

        //- Point of hit; invalid for misses
        Point hitPoint_;

        //- Label of face hit
        label index_;


public:

    // Constructors

        //- Construct from components
        PointIndexHit(const bool success, const Point& p, const label index)
        :
            hit_(success),
            hitPoint_(p),
            index_(index)
        {}

        //- Construct from point. Hit and distance set later
        PointIndexHit(const Point& p)
        :
            hit_(false),
            hitPoint_(p),
            index_(-1)
        {}

        //- Construct null
        PointIndexHit()
        :
            hit_(false),
            hitPoint_(Zero),
            index_(-1)
        {}

        //- Construct from Istream
        PointIndexHit(Istream& is)
        {
            is >> *this;
        }


    // Member Functions

        //- Is there a hit
        bool hit() const
        {
            return hit_;
        }

        //- Return index
        label index() const
        {
            return index_;
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

        Point& rawPoint()
        {
            return hitPoint_;
        }

        void setHit()
        {
            hit_ = true;
        }

        void setMiss()
        {
            hit_ = false;
        }

        void setPoint(const Point& p)
        {
            hitPoint_ = p;
        }

        void setIndex(const label index)
        {
            index_ = index;
        }

        bool operator==(const PointIndexHit& rhs) const
        {
            return
                hit_ == rhs.hit()
             && hitPoint_ == rhs.rawPoint()
             && index_ == rhs.index();
        }

        bool operator!=(const PointIndexHit& rhs) const
        {
            return !operator==(rhs);
        }

        void write(Ostream& os)
        {
            if (hit())
            {
                os << "hit:" << hitPoint() << " index:" << index();
            }
            else
            {
                os << "miss:" << missPoint() << " index:" << index();
            }
        }

        friend Ostream& operator<< (Ostream& os, const PointIndexHit& pHit)
        {
            if (os.format() == IOstream::ASCII)
            {
                os  << pHit.hit_ << token::SPACE << pHit.hitPoint_
                    << token::SPACE << pHit.index_;
            }
            else
            {
                os.write
                (
                    reinterpret_cast<const char*>(&pHit),
                    sizeof(PointIndexHit)
                );
            }

            // Check state of Ostream
            os.check("Ostream& operator<<(Ostream&, const PointIndexHit&)");

            return os;
        }

        friend Istream& operator>>(Istream& is, PointIndexHit& pHit)
        {
            if (is.format() == IOstream::ASCII)
            {
                return is >> pHit.hit_ >> pHit.hitPoint_ >> pHit.index_;
            }
            else
            {
                is.read
                (
                    reinterpret_cast<char*>(&pHit),
                    sizeof(PointIndexHit)
                );
            }

            // Check state of Istream
            is.check("Istream& operator>>(Istream&, PointIndexHit&)");

            return is;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
