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

Class
    Foam::pTraits

Description
    Traits class for primitives.

    All primitives need a specialised version of this class. The
    specialised version will normally also require a conversion
    method.

\*---------------------------------------------------------------------------*/

#ifndef pTraits_H
#define pTraits_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;

/*---------------------------------------------------------------------------*\
                           Class pTraits Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class pTraits
:
    public PrimitiveType
{

public:

    // Constructors

        //- Construct from primitive
        pTraits(const PrimitiveType& p)
        :
            PrimitiveType(p)
        {}

        //- Construct from Istream
        pTraits(Istream& is)
        :
            PrimitiveType(is)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
