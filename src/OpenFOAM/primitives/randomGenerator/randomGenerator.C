/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "randomGenerator.H"
#include "uint64.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomGenerator::randomGenerator(Istream& is)
:
    x_(pTraits<uint64_t>(is))
{}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, randomGenerator& rndGen)
{
    is >> rndGen.x_;
    is.check("operator>>(Istream& is, randomGenerator& rndGen)");
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const randomGenerator& rndGen)
{
    os << rndGen.x_;
    os.check("operator<<(Ostream& os, const randomGenerator& rndGen)");
    return os;
}


// ************************************************************************* //
