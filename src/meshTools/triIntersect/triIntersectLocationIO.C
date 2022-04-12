/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "triIntersectLocation.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * IOstream Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::triIntersect::operator<<
(
    Ostream& os,
    const location& l
)
{
    if (l.isSrcPoint())
    {
        os << "s" << l.srcPointi();
    }
    if (l.isSrcPoint() && l.isTgtPoint())
    {
        os << '&';
    }
    if (l.isTgtPoint())
    {
        os << "t" << l.tgtPointi();
    }
    if (l.isIntersection())
    {
        os << "s" << l.srcEdgei() << "^t" << l.tgtEdgei();
    }
    return os;
}


// ************************************************************************* //
