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
    Foam::volumeType

Description

SourceFiles
    volumeType.C

\*---------------------------------------------------------------------------*/

#ifndef volumeType_H
#define volumeType_H

#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class volumeType;
Istream& operator>>(Istream& is, volumeType&);
Ostream& operator<<(Ostream& os, const volumeType& C);


/*---------------------------------------------------------------------------*\
                         Class volumeType Declaration
\*---------------------------------------------------------------------------*/

class volumeType
{
public:

    //- Volume types
    enum type
    {
        UNKNOWN = 0,
        MIXED   = 1,
        INSIDE  = 2,
        OUTSIDE = 3
    };


private:

    // Private data

        //- Volume type
        type t_;


public:

    // Static data

        static const NamedEnum<volumeType, 4> names;


    // Constructors

        //- Construct null
        volumeType()
        :
            t_(UNKNOWN)
        {}

        //- Construct from components
        volumeType(type t)
        :
            t_(t)
        {}


    // Member Functions

        operator type() const
        {
            return t_;
        }


    // IOstream operators

        friend Istream& operator>>(Istream& is, volumeType& vt);
        friend Ostream& operator<<(Ostream& os, const volumeType& vt);
};


//- Data associated with volumeType type are contiguous
template<>
inline bool contiguous<volumeType>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
