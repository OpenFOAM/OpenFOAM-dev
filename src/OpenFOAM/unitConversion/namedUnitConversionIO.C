/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "namedUnitConversion.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::namedUnitConversion::readName(const entry& e)
{
    // Read the entry characters directly
    OStringStream oss;
    dynamicCast<const primitiveEntry>(e).write(oss, true);
    name_ = oss.str();

    // Beautify slightly by removing unnecessary spaces
    static const string a(List<char>({token::BEGIN_SQR, token::SPACE}));
    static const string b(List<char>({token::BEGIN_SQR}));
    name_.replace(a, b);
}


void Foam::namedUnitConversion::readName
(
    const word& keyword,
    const dictionary& dict
)
{
    readName(dict.lookupEntry(keyword, false, true));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::namedUnitConversion::read
(
    const word& keyword,
    const dictionary& dict
)
{
    unitConversion::read(keyword, dict);

    readName(keyword, dict);
}


bool Foam::namedUnitConversion::readIfPresent
(
    const word& keyword,
    const dictionary& dict
)
{
    if (unitConversion::readIfPresent(keyword, dict))
    {
        readName(keyword, dict);
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, namedUnitConversion& u)
{
    const primitiveEntry e(word::null, is);
    u.reset(unitConversion(e.stream()));
    u.readName(e);

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, namedUnitConversion&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const namedUnitConversion& u)
{
    os << u.name_.c_str();

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const namedUnitConversion&)");

    return os;
}


// ************************************************************************* //
