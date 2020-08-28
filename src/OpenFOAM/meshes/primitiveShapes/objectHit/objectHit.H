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
    Foam::objectHit

Description
    This class describes a combination of target object index and success flag.

\*---------------------------------------------------------------------------*/

#ifndef objectHit_H
#define objectHit_H

#include "bool.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class objectHit Declaration
\*---------------------------------------------------------------------------*/

class objectHit
{
    // Private Data

        //- Hit success
        bool hit_;

        //- Object of hit
        label hitObject_;


public:

    // Constructors

        //- Construct null
        objectHit()
        :
            hit_(false),
            hitObject_(-1)
        {}

        //- Construct from components
        objectHit(const bool success, const label obj)
        :
            hit_(success),
            hitObject_(obj)
        {}

        //- Construct from Istream
        objectHit(Istream& is)
        :
            hit_(readBool(is)),
            hitObject_(readLabel(is))
        {}


    // Member Functions

        //- Is there a hit
        bool hit() const
        {
            return hit_;
        }

        //- Return hit object
        label hitObject() const
        {
            return hitObject_;
        }


    // Friend Operators

        inline friend bool operator==(const objectHit& a, const objectHit& b)
        {
            return (a.hit_ == b.hit_) && (a.hitObject_ == b.hitObject_);
        }

        inline friend bool operator!=(const objectHit& a, const objectHit& b)
        {
            return !(a == b);
        }


    // Ostream operator

        inline friend Ostream& operator<<(Ostream& os, const objectHit& obj)
        {
             return os << obj.hit() << token::SPACE << obj.hitObject();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
