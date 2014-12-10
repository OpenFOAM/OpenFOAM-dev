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

Class
    Foam::labelBits

Description
    A 29bits label and 3bits direction packed into single label

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef labelBits_H
#define labelBits_H

#include "label.H"
//#include "uLabel.H"
#include "direction.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class labelBits Declaration
\*---------------------------------------------------------------------------*/

class labelBits
{
    // Private data

        label data_;

        inline static label pack(const label val, const direction bits)
        {
#           ifdef FULLDEBUG
            if (bits > 7 || (((val<<3)>>3) != val))
            {
                FatalErrorIn
                (
                    "labelBits::pack(const label, const direction)"
                )   << "Direction " << bits << " outside range 0..7"
                    << " or value " << val << " negative or larger than "
                    << label(8*sizeof(label)-3) << " bit representation"
                    << abort(FatalError);
            }
#           endif

            return (val<<3) | bits;
        }

public:

    // Constructors

        //- Construct null
        inline labelBits()
        {}

        //- Construct from components
        inline labelBits(const label val, const direction bits)
        :
             data_(pack(val, bits))
        {}

        //- Construct from Istream
        inline labelBits(Istream& is)
        {
            is >> data_;
        }



    // Member Functions

        inline label val() const
        {
            return data_ >> 3;
        }

        inline direction bits() const
        {
            return data_ & 7;
        }

        inline void setVal(const label val)
        {
            data_ = pack(val, bits());
        }

        inline void setBits(const direction bits)
        {
            data_ = pack(val(), bits);
        }


    // Member Operators

        inline friend bool operator==(const labelBits& a, const labelBits& b)
        {
            return a.data_ == b.data_;
        }

        inline friend bool operator!=(const labelBits& a, const labelBits& b)
        {
            return !(a == b);
        }

    // IOstream Operators

        inline friend Istream& operator>>(Istream& is, labelBits& lb)
        {
            return is >> lb.data_;
        }

        inline friend Ostream& operator<<(Ostream& os, const labelBits& lb)
        {
            return os << lb.data_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
