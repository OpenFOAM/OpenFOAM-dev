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

#include "pointData.H"

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointData& wDist)
{
    if (os.format() == IOstream::ASCII)
    {
        return os
            << static_cast<const pointEdgePoint&>(wDist)
            << token::SPACE << wDist.s() << token::SPACE << wDist.v();
    }
    else
    {
        return os
            << static_cast<const pointEdgePoint&>(wDist)
            << wDist.s() << wDist.v();
    }
}


Foam::Istream& Foam::operator>>(Istream& is, pointData& wDist)
{
    return is >> static_cast<pointEdgePoint&>(wDist) >> wDist.s_ >> wDist.v_;
}


// ************************************************************************* //
