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

Foam::randomGenerator::randomGenerator(Istream& is, const bool global)
:
    global_(global),
    x_(pTraits<uint64_t>(is))
{
    checkSync();
}


Foam::randomGenerator::randomGenerator
(
    const word& name,
    const dictionary& dict,
    randomGenerator&& defaultRndGen
)
:
    global_(defaultRndGen.global_),
    x_
    (
        dict.found(name)
      ? dict.lookup<uint64_t>(name)
      : dict.found(name + "Seed")
      ? seed(dict.lookup<label>(name + "Seed")).x(global_)
      : defaultRndGen.x_
    )
{}


Foam::randomGenerator::randomGenerator
(
    const word& name,
    const dictionary& dict,
    const seed defaultS,
    const bool global
)
:
    randomGenerator(name, dict, randomGenerator(defaultS, global))
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::randomGenerator::operator=(const randomGenerator& rndGen)
{
    if (global_ != rndGen.global_)
    {
        FatalErrorInFunction
            << "Attempted assignment of a " << (global_ ? "" : "non-")
            << "global random generator to a " << (rndGen.global_ ? "" : "non-")
            << "global random generator"
            << exit(FatalError);
    }
    x_ = rndGen.x_;
    checkSync();
}


void Foam::randomGenerator::operator=(randomGenerator&& rndGen)
{
    *this = static_cast<const randomGenerator&>(rndGen);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, randomGenerator& rndGen)
{
    is >> rndGen.x_;
    rndGen.checkSync();
    is.check("operator>>(Istream& is, randomGenerator& rndGen)");
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const randomGenerator& rndGen)
{
    rndGen.checkSync();
    os << rndGen.x_;
    os.check("operator<<(Ostream& os, const randomGenerator& rndGen)");
    return os;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::writeEntry(Ostream& os, const randomGenerator& rndGen)
{
    os << rndGen;
}


// ************************************************************************* //
